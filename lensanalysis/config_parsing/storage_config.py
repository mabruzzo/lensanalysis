import ConfigParser

from ..misc.serialization import ShearMapCollectionFGStorage, \
    ConvergenceMapCollectionFGStorage, PeakLocCollectionFGStorage, \
    PeakCountCollectionFGStorage, FullShearCatFGLoader

from ..misc.fname_formatter import BaseFnameFormatter, \
    RealizationBinnedSubdirectoryFormatter

def _set_defaults_dynamically(config, section, defaults):
    for elem in defaults:
        try:
            config.get(section,elem[0])
        except ConfigParser.NoOptionError:
            if isinstance(elem[1],str):
                val = elem[1]
            else:
                val = str(elem[1])
            config.set(section,elem[0],val)

_tomo_descriptors = ["tomo","tomography","tomographic"]
_smoothing_descriptors = ['smoothed','smooth']
_noisy_descriptors = ["noisy"]
            
def _check_unallowed_descriptors(descriptors,descriptor_subsets):
    """
    Helper function that ensures descriptors only fall in a subset.
    """
    for elem in descriptor_subsets:
        if elem == 'noisy':
            ref_descriptors = _noisy_descriptors
        elif elem == 'tomo':
            ref_descriptors = _tomo_descriptors
        elif elem == 'smooth':
            ref_descriptors = _smoothing_descriptors
        else:
            raise ValueError(("{:s} is invalid descriptor checking"
                              "subset").format(elem))
        for i in ref_descriptors:
            if i in descriptors:
                raise ValueError("{:s} is not an allowed descriptor".format(i))
                
        
def _parse_descriptor(descriptors):
    """
    Helper Function.
    """
    if len(descriptors) > 3:
        raise ValueError("Too many descriptors mentioned")
    results = [False,False, False]
    for elem in descriptors:
        if elem.lower() in _tomo_descriptors:
            if results[0]:
                raise ValueError("tomographic descriptor mentioned more than "
                                 "once.")
            results[0] = True
            raise NotImplementedError("Not Presently Equiped to handle "
                                      "tomographic data")

        elif elem.lower() in _smoothing_descriptors:
            if results[1]:
                raise ValueError("smoothing descriptor mentioned more than "
                                 "once.")
            results[1] = True

        elif elem.lower() in _noisy_descriptors:
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
        temp_l.append('noisy')
    return '_'.join(temp_l)

def _check_create_storage(config,section,option,default = None):
    """
    Checks to see if the configuration file want to ever allow this storage.

    Parameters
    ----------
    config - ConfigParser
        Instance of RawConfigParser, ConfigParser or SafeConfigParser
    section - str
        The name of the section where the configuration comes from
    option - str
        The name of the option that holds this info
    default - None or bool, optional
        The default value to return if there is no such option. If the value of 
        this argument is None, then an error is raised if the option does not 
        exist. By default this is set to None
    """

    if not config.has_section(section):
        message = "The [{:s}] section is missing.".format(section)
        raise ConfigParser.NoSectionError(message)

    if not config.has_option(section,option):
        if default is not None:
            # if default is None we will allow for the error
            return default
    return config.getboolean(section,option)

def _option_builder(*args):
    temp = [arg for arg in args if arg is not None]
    return '_'.join(temp)

def _create_collection_storage(descriptors,root_dir,section_name,
                               config, core_option_name = None,
                               storage_option_suffix = None,
                               storage_dir_suffix = 'dir',
                               storage_fname_suffix = 'fname',
                               noiseless_prefix = True):
    """
    General Helper Function to help read in storage configuration.

    Parameters
    ----------
    descriptors : sequence or str
        descriptors to discriminate between types of storage
    root_dir : str
        the directory within which all data is stored
    section - str
        The name of the section where the configuration comes from
    config : ConfigParser
        Instance of RawConfigParser, ConfigParser or SafeConfigParser
    core_option_name : str or None, optional
        The core part of the option name. If this is None, then there is no 
        core part of the string. By default this is None.
    storage_option_suffix : str or None
        The suffix of the option that tells you whether or not we want to make 
        a storage object. By default this is None.
    storage_dir_suffix : str or None
        The suffix of the option that tells you where the name of the storage 
        directory that will (or already does) hold the collections. Default is 
        'dir'.
    storage_fname_suffix : str or None
        The suffix of the option that gives the file name template. Default is 
        'fname'.
    """

    decoded_arg = _parse_descriptor(descriptors)
    string_prefix = _descriptor_string_prefix(decoded_arg)

    if not noiseless_prefix and string_prefix == "noiseless":
        string_prefix = None
    make_storage_option = _option_builder(string_prefix, core_option_name,
                                          storage_option_suffix)
    make_storage = _check_create_storage(config, section_name,
                                         make_storage_option,
                                         default = False)
    if not make_storage:
        return None

    # Now that we know that storage is allowed, check the remaining fields

    # find the local directory where the storage is held
    dir_option_name = _option_builder(string_prefix, core_option_name,
                                      storage_dir_suffix)
    local_dir = config.get(section_name, dir_option_name)
    storage_dir = '/'.join((root_dir,local_dir))

    # Find the fname_template
    fname_option_name = _option_builder(string_prefix, core_option_name,
                                        storage_fname_suffix)
    fname_template = config.get(section_name, fname_option_name)

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

    We should almost certainly be copying over the data from the ConfigParser 
    and ignoring the ConfigParser afterwards.
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
                                          self._config_parser,
                                          core_option_name = None,
                                          storage_option_suffix = 'map',
                                          storage_dir_suffix = 'dir',
                                          storage_fname_suffix = 'fname')

    def shear_map_collection_storage(self,descriptors,root_dir):
        _check_unallowed_descriptors(descriptors,['smooth'])
        return _create_collection_storage(descriptors, root_dir,
                                          "ShearMaps",
                                          self._config_parser,
                                          core_option_name = None,
                                          storage_option_suffix = 'map',
                                          storage_dir_suffix = 'dir',
                                          storage_fname_suffix = 'fname')

    def peak_loc_collection_storage(self,descriptors,root_dir):
        _check_unallowed_descriptors(descriptors,['smooth','noisy'])
        return _create_collection_storage(descriptors, root_dir,
                                          "FeatureProducts",
                                          self._config_parser,
                                          core_option_name = "peak_loc",
                                          storage_option_suffix = None,
                                          storage_dir_suffix = 'dir',
                                          storage_fname_suffix = 'fname',
                                          noiseless_prefix = False)

    def peak_counts_collection_storage(self,descriptors,root_dir):
        _check_unallowed_descriptors(descriptors,['smooth','noisy'])
        return _create_collection_storage(descriptors, root_dir,
                                          "FeatureProducts",
                                          self._config_parser,
                                          core_option_name = "peak_counts",
                                          storage_option_suffix = None,
                                          storage_dir_suffix = 'dir',
                                          storage_fname_suffix = 'fname',
                                          noiseless_prefix = False)


def _set_subdir_binned_realization_defaults(config, section, prefix):
    """
    We set the defaults in this way rather than the built-in way because to use 
    the built-in method we most know the prefixes ahead of time.
    """
    defaults = [('_'.join((prefix,'subdir')), False),
                ('_'.join((prefix,'subdir_min_loc')), 
                 'binned_realizations'),
                ('_'.join((prefix,'subdir_min_loc')), -1),
                ('_'.join((prefix,'subdir_max_loc')), -1),
                ('_'.join((prefix,'subdir_mid_loc')), -1)]
    _set_defaults_dynamically(config,section,defaults)

def _build_subdir_binned_realization_formatter(fname_formatter, config, section,
                                               prefix):
    _set_subdir_binned_realization_defaults(config,section,prefix)
    subdir_choice = '_'.join((prefix,'subdir'))

    if not config.getboolean(section,subdir_choice):
        return fname_formatter

    format_option = config.get(section,'_'.join((prefix,'subdir_format')))
    assert format_option == 'binned_realizations'
    
    directory_template = config.get(section,'_'.join((prefix,
                                                      'subdir_template')))
    realizations_per_bin=config.getint(section,
                                       '_'.join((prefix,
                                                 'subdir_realizations_per_bin')))

    fields = [None, None, None]

    iterable = [(config.getint(section,'_'.join((prefix,'subdir_min_loc'))),
                 'min'),
                (config.getint(section,'_'.join((prefix,'subdir_max_loc'))),
                 'max'),
                (config.getint(section,'_'.join((prefix,'subdir_mid_loc'))),
                 'mid')]
    for field_loc, field in iterable:
        
        if field_loc != -1:
            if 0<=field_loc < 3:
                if fields[field_loc] is None:
                    fields[field_loc] = 'min'
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
                                                  realizations_per_bin,
                                                  fname_formatter,
                                                  "realization")

def _construct_shear_fname_formatter(config, section):
    fname_template = config.get(section,"shear_cat_fname_template")
    bin_loc = config.getint(section,"shear_cat_fname_binning_loc")
    if bin_loc not in [0,1]:
        raise ValueError("shear_cat_fname_binning_loc must be 0 or 1")
    realization_loc = config.getint(section,
                                    "shear_cat_fname_realization_loc")
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
        pos_formatter = _pos_fname_formatter(self._config,"ShearCats")
        loader = FullShearCatFGLoader(shear_formatter, root_dir, num_elements,
                                      _normal_field_mapping, pos_formatter)
        return loader

if __name__ == '__main__':
    ShearCatCollectionLoaderConfig("shear_cat_config.ini")
