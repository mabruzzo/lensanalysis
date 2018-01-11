import ConfigParser

from ..misc.serialization import ShearMapCollectionFGStorage, \
    ConvergenceMapCollectionFGStorage, PeakLocCollectionFGStorage, \
    PeakCountCollectionFGStorage, BaseFnameFormatter

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
    # check to see if the storage directory exists
    raise RuntimeError("CHECK IF STORAGE DIRECTORY EXISTS - if it does "
                       "create it (within reason).")

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
        

