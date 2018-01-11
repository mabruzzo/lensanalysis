import ConfigParser

from ..misc.serialization import ShearMapCollectionFGStorage, \
    ConvergenceMapCollectionFGStorage, PeakLocCollectionFGStorage, \
    PeakCountCollectionFGStorage, FullShearCatFGLoader, BaseFnameFormatter, \
    RealizationBinnedSubdirectoryFormatter

def _parse_descriptor(descriptors):
    """
    Helper Function.
    """
    if len(descriptors) > 3:
        raise ValueError("Too many descriptors mentioned")
    results = [False,False, False]
    for elem in descriptors:
        if elem.lower() in ["tomo","tomography","tomographic"]:
            if results[0]:
                raise ValueError("tomographic descriptor mentioned more than "
                                 "once.")
            results[0] = True
            raise NotImplementedError("Not Presently Equiped to handle "
                                      "tomographic data")

        elif elem.lower() in ['smoothed','smooth']:
            if results[1]:
                raise ValueError("smoothing descriptor mentioned more than "
                                 "once.")
            results[1] = True

        elif elem.lower() == "noisy":
            if results[2]:
                raise ValueError("noisy descriptor mentioned more than "
                                 "once.")
            results[2] = True

        elif elem.lower() != 'noiseless':
            raise ValueError("{:s} is not an allowed descriptor.".format(elem))

    return results

def _descriptor_string_prefix(results):
    if not any(results):
        return "noiseless"
    elif results[0] and not any(results[1:]):
        return "tomo_noiseless"

    temp_l = []

    if results[0]:
        temp_l.append('tomo')
    if results[1]:
        temp_l.append('smoothed')
    if results[2]:
        temp_l.append('noisy'):
    return '_'.join(temp_l)

def _create_collection_storage(descriptors,root_dir,section_name,config):
    decoded_arg = _parse_descriptor(descriptors)
    string_prefix = _descriptor_string_prefix

    config= self._config_parser
        
    # check to see if ConvergenceMaps was a section in the configuration
    # if not then we won't save anything.
    if not config.has_section(section_name):
        return None

    # probably should set up some formal defaults

    # check to see if the line allowing for storage creation section exists
    # (if it does not assume that it is not allowed)
    if not config.has_section(section_name, '_'.join([string_prefix,'map'])):
            return None

    if config.has_option(section_name,'_'.join([string_prefix,'map'])):
        return None
    make_storage = config.getboolean(section_name,
                                     '_'.join([string_prefix,'map']))
    if not make_storage:
        return None

    # Now that we know that storage is allowed, check the remaining fields

    # find the local directory where the storage is held
    local_dir = config.get(section_name, '_'.join([string_prefix,'dir']))
    storage_dir = '/'.join(root_dir,local_dir)

    # Find the fname_template
    fname_template = config.get(section_name,
                                '_'.join([string_prefix,'fname']))

    field_mapping = {"collection_id":"realization"}

    if decoded_arg[0]:
        raise NotImplementedError("Not currently allowing for the creation "
                                  "of tomographic storage. We need to "
                                  "create the machinary for the "
                                  "subdirectories.")
    else:
        fname_formatter = BaseFnameFormatter(fname_template,"realization")
        result = ConvergenceMapCollectionFGStorage(fname_formatter,
                                                   storage_dir, 1,
                                                   field_mapping)
    return result

class StorageConfig(object):
    """
    Represents the storage configuration specified by a configuration file.
    """

    def __init__(self,config_file):
        self._config_parser = ConfigParser.SafeConfigParser()
        self._config_parser.read(config_file)

    def convergence_map_collection_storage(self,descriptors,root_dir):
        """
        Returns instance of ConvergenceStorage if specified in the storage 
        configuration. If not specified, then None is returned.
        """
        return _create_collection_storage(descriptors, root_dir,
                                          "ConvergenceMaps",
                                          self._config_parser)

    def shear_map_collection_storage(self,descriptors,root_dir):
        if "smooth" in descriptors:
            raise ValueError("smooth is not an allowed descriptor")
        if "smoothed" in descriptors:
            raise ValueError("smoothed is not an allowed descriptor")
        return _create_collection_storage(descriptors, root_dir,
                                          "ShearMaps", self._config_parser)

    def peak_loc_collection_storage(self,descriptors,root_dir):
        pass

    def peak_counts_collection_storage(self,descriptors,root_dir):
        pass
        

def _set_defaults_dynamically(config, section, defaults):
    for elem in defaults:
        try:
            config.get(section,elem[0])
        except ConfigParser.NoOptionError:
            config.set(section,elem[0],elem[1])
    
def _set_subdir_binned_realization_defaults(config, section, prefix):
    """
    We set the defaults in this way rather than the built-in way because to use 
    the built-in method we most know the prefixes ahead of time.
    """
    defaults = [('_'.join(prefix,'subdir'), False),
                ('_'.join(prefix,'subdir_min_loc'), 'binned_realizations'),
                ('_'.join(prefix,'subdir_min_loc'), -1),
                ('_'.join(prefix,'subdir_max_loc'), -1),
                ('_'.join(prefix,'subdir_mid_loc'), -1)]
    _set_defaults_dynamically(config,section,defaults)

def _build_subdir_binned_realization_formatter(fname_formatter, config, section,
                                               prefix):
    _set_subdir_binned_realization_defaults(config,section,prefix)
    subdir_choice = '_'.join(prefix,'subdir')

    if not config.getboolean(section,elem[0]):
        return fname_formatter

    format_option = config.get(section,'_'.join(prefix,'subdir_min_loc')) 
    assert format_option == 'binned_realizations'
    
    directory_template = config.get(section,'_'.join(prefix,'subdir_template'))
    realizations_per_bin=config.getint(section,
                                       '_'.join(prefix,
                                                'subdir_realizations_per_bin'))

    fields = [None, None, None]

    iterable = [(config.getint(section,'_'.join(prefix,'subdir_min_loc')),
                 'min'),
                (config.getint(section,'_'.join(prefix,'subdir_max_loc')),
                 'max'),
                (config.getint(section,'_'.join(prefix,'subdir_mid_loc')),
                 'mid')]
    for field_loc, field in iterable:
        
        if min_loc != -1:
            if 0<=min_loc < 3:
                if fields[min_loc] is None:
                    fields[min_loc] = 'min'
                else:
                    raise ValueError("More than one field specified for the "
                                     "same location in the template")
            else:
                raise ValueError(("{:s} must be assigned an integer value "
                                  "location from -1 to 2").format(field))
    if fields[2] is None:
        if fields[1] is None:
            if fields[0] is None:
                raise ValueError("A subdirectory is unecessary")
            else:
                fields = fields[:1]
        else:
            fields = fields[:2]
    return RealizationBinnedSubdirectoryFormatter(directory_template,
                                                  fields,
                                                  realizations_per_dir,
                                                  fname_formatter,
                                                  "realization")

def _construct_shear_fname_formatter(config, section):
    fname_template = config.get(section,"shear_cat_fname_template")
    bin_loc = config.get(section,"shear_cat_fname_binning_loc")
    if bin_loc not in [0,1]:
        raise ValueError("shear_cat_fname_binning_loc must be 0 or 1")
    realization_loc = config.get(section,"shear_cat_fname_realization_loc")
    if realization_loc not in [0,1]:
        raise ValueError("shear_cat_fname_realization_loc must be 0 or 1")
    if realization_loc == bin_loc:
        raise ValueError("shear_cat_fname_realization_loc and "
                         "shear_cat_fname_binning_loc must\nhave different "
                         "values")
    if bin_loc == 0:
        fields = ["bin", "realization"]
    else:
        fields = ["realization", "bin"]
    fname_formatter = BaseFnameFormatter(fname_template, fields)
    return _build_subdir_binned_realization_formatter(fname_formatter, config,
                                                      section, 'shear_cat')
def _pos_fname_formatter(config,section):
    fname_template = config.get(section,"shear_cat_fname_template")
    return BaseFnameFormatter(fname_template, ["bin"])

_normal_field_mapping = {"collection_id" : "realization",
                         "element_id" : "bin"}

class ShearCatCollectionLoaderConfig(object):
    """
    Object responsible for creating the shear catalog loader.
    """

    def __init__(self,config_file):
        config = ConfigParser.SafeConfigParser()
        config.read(config_file)
        self._config = config

    def constructCompleteShearCatLoader(self,root_dir):
        num_elements = self._config.getint("ShearCats","num_cat_bins")
        shear_formatter = _construct_shear_fname_formatter(self._config,
                                                           "ShearCats")
        pos_formatter = _pos_fname_formatter(config,"ShearCats")
        loader = FullShearCatFGLoader(shear_formatter, root_dir, num_elements,
                                      _normal_field_mapping, pos_formatter)
        return loader

if __name__ == '__main__':
    ShearCatCollectionLoaderConfig("shear_cat_config.ini")
