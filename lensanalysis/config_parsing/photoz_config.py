import ConfigParser
import re

from ..procedure.tomo_binning import ConstantBias, InvertedConstantBias

def _option_reader(config,section_name,required_option_regex):
    options = config.options(section_name)

    length = len(options)
    i = 0

    option_nums = []
    option_names = []

    while i < len(options):
        option = options[i]
        match = required_option_regex.match(option)
        if match is None:
            i+=1
            continue
        integer = int(match.groups()[0])
        if integer not in option_nums:
            option_nums.append(integer)
            option_names.append(option)
        else:
            raise ValueError("{:s} was specified multiple times".format(option))
        options.remove(option)
    return options,option_nums,option_names

def _retrieval_method(config,type=None):
    if type is None:
        method = config.get
    elif type == float:
        method = config.getfloat
    elif type == int:
        method = config.getint
    elif type == bool:
        method = config.getboolean
    else:
        raise ValueError("type must be None, float or int")
    return method

def _build_map(config,section_name,required_option_regex,type=None):
    option_nums,option_names = _option_reader(config,section_name,
                                              required_option_regex)[1:]

    out = {}
    method = _retrieval_method(config,type)
    
    for option_num,option_name in zip(option_nums,option_names):
        option_val = method(section_name,option_name)
        if option_val in out:
            raise ValueError("More than one option is identified as "
                             "{:s}".format(option_val))
        out[option_val] = option_num
    return out

class VariableSectionReader(object):
    """
    option_re is a RegexObject which has a single group. That group must be 
    convertible to an integer.
    """
    def __init__(self,config,section_name,required_option_regex,
                 other_options=[],forwarding_section = None,
                 identifier_type=None):
        self.config = config
        self.section_name = section_name
        self.required_option_regex = required_option_regex
        self.other_options = other_options
        self.forwarding_section = forwarding_section

        self._identity_map = _build_map(config,section_name,
                                        required_option_regex,
                                        type=identifier_type)

    def get_option_nums(self):
        return self._identity_map.values()

    def get_identifiers(self):
        return self._identity_map.keys()

    def get_identifier_num(self,identifier):
        return self._identity_map[identifier]

    def check_correctness(self):
        options,option_nums = _option_reader(self.config,self.section_name,
                                             self.required_option_regex)[:2]

        for option_num in option_nums:
            for option_name, option_template,val_type in self.other_options:
                option_name = elem.format(option_num)
                try:
                    options.remove(option_name)
                except ValueError:
                    pass
        return len(options) == 0

    def build_info_dict(self,identifier):
        num = self.get_identifier_num(identifier)

        out = {}
        for option_name, option_template,val_type in self.other_options:
            method = _retrieval_method(self.config,type=val_type)
            
            try:
                val = method(self.section_name,option_template.format(num))
            except ConfigParser.NoOptionError:
                if self.forwarding_section is not None:
                    val = method(self.forwarding_section,option_name)
                else:
                    continue
            out[option_name] = val
        return out

    def build_all_info_dicts(self,config):
        out = {}
        for identifier in self._identity_map.iterkeys():
            out[identifier] = self.build_info_dict(identifier)
        return out


#_photo_z_section_reader = VariableSectionReader(
_photo_z_options = [("bias_factor","bias_factor_{:d}",float),
                    ("scatter_factor","scatter_factor_{:d}",float),
                    ("stochastic_seed","stochastic_seed_{:d}",int),
                    ("constant_seed","constant_seed_{:d}",bool),
                    ("zphot_start", "zphot_start_{:d}", bool)]

class PhotozConfig(object):
    def __init__(self,config):
        self._config=config
        self._section_reader = VariableSectionReader(config,"PhotozDataset",
                                                     re.compile('name_(.*)'),
                                                     _photo_z_options,
                                                     "StandardSettings",
                                                     None)

    def has_identifier(self,identifier):
        return identifier in self._section_reader.get_identifiers()

    def get_photoz_noise_addition(self,identifier):
        raise NotImplementedError()

class PseudoPhotozConfig(PhotozConfig):
    def has_identifier(self,identifier):
        # remove the "_ppz" suffix
        return identifier[:-4] in self._section_reader.get_identifiers()

    def get_fid_name(self,identifier):
        num = self._section_reader.get_identifier_num(identifier[:-4])
        try:
            return self._config.get("PseudoPhotozSource",
                                    "fid_name_{:d}".format(num))
        except ConfigParser.NoOptionError:
            return None

    def get_photoz_noise_addition(self,identifier):
        photoz_info = self._section_reader.build_info_dict(identifier[:-4])

        # for now we only allow constant bias.
        if photoz_info["scatter_factor"] == 0:
            if photoz_info["zphot_start"]:
                return InvertedConstantBias(photoz_info["bias_factor"])
            else:
                return ConstantBias(photoz_info["bias_factor"])
        else:
            raise NotImplementedError("Not yet equipped to handle stochastic "
                                      "noise")
