import ConfigParser
from collections import OrderedDict
import os
import os.path

from .storage_config import StorageConfig, ShearCatCollectionLoaderConfig

from ..misc.analysis_collection import AnalysisProductCollection, \
    ConvergenceMapProductCollection, ShearMapProductCollection, \
    FeatureProductCollection

def _load_single_storage_collection(path, storage_config_method, descriptions):
    """
    This is probably not the way to do this.
    """
    storage = storage_config_method(descriptions,path)
    if storage is None:
        return None
    new_path = os.path.normpath(storage.root_dir)
    if not os.path.isdir(new_path):
        # let's try to make a new path
        head,tail = os.path.split(new_path)
        if os.path.isdir(head):
            os.mkdir(new_path)
            if not os.path.isdir(new_path):
                raise RuntimeError("IT SHOULDN'T BE POSSIBLE FOR THIS TO GET "
                                   "CALLED")
        else:
            raise ValueError("{:s} does not exist/ is not a directory")
    return storage

def _load_full_storage_collection(path,storage_config):
    """
    Create the storage collection as you go.

    For now let's just build everything we can. Probably not a great idea.

    Come back and update after allowing for tomography.
    """
    analysis_collection = AnalysisProductCollection()

    # first handle non-tomographic convergence maps
    # kappa storage collection
    ksc = ConvergenceMapProductCollection()
    storage_method = storage_config.convergence_map_collection_storage
    ksc.noiseless_map = _load_single_storage_collection(path,
                                                        storage_method,
                                                        [])
    ksc.noisy_map = _load_single_storage_collection(path,
                                                    storage_method,
                                                    ['noisy'])
    ksc.smoothed_map = _load_single_storage_collection(path,
                                                       storage_method,
                                                       ['smoothed'])
    ksc.smoothed_noisy_map = _load_single_storage_collection(path,
                                                             storage_method,
                                                             ['smoothed',
                                                              'noisy'])
    analysis_collection.conv_map=ksc

    # shear storage collection
    ssc = ShearMapProductCollection()
    storage_method = storage_config.shear_map_collection_storage
    
    ssc.noiseless_map = _load_single_storage_collection(path,
                                                        storage_method,
                                                        [])
    ssc.noisy_map = _load_single_storage_collection(path,
                                                    storage_method,
                                                    ['noisy'])
    analysis_collection.shear_map = ssc

    # now, finally let's load in the feature product collection
    fpc = FeatureProductCollection()

    storage_method = storage_config.peak_loc_collection_storage
    fpc.peak_locations = _load_single_storage_collection(path,
                                                         storage_method,
                                                         [])
    storage_method = storage_config.peak_counts_collection_storage
    fpc.peak_counts = _load_single_storage_collection(path,
                                                      storage_method,
                                                      [])
    analysis_collection.feature_products=fpc
    return analysis_collection


class Cache(object):
    def __init__(self):
        self._dict = {}

    def __getitem__(self,key):
        return self._dict[key]
        
    def __setitem__(self,key,item):
        self._dict[key] = item

    def keys(self):
        return self._dict.keys()

    def __contains__(self,item):
        return item in self._dict

class SizeLimitedCache(Cache):
    def __init__(self,max_items):
        self._max_items = max_items
        self._dict = OrderedDict()

    def __setitem__(self,key,item):
        if len(self._dict) == self._maxitems:
            self._dict.popitem(last=True)
        self._dict[key] = item

    def __getitem__(self,key):
        val = self._dict[key]
        # update access order
        del self._dict[key]
        self._dict[key] = val
        return val


class CosmologyAnalysisCollection(object):
    """
    For now we will just base the collection off of the directory names.

    This class probably has too much responsibility.
    """

    def __init__(self,root, storage_config, root_shear, shear_cat_config,
                 cache = None):
        root_path = os.path.normpath(root)
        if not os.path.exists(root_path):
            raise ValueError("{:s} does not already exist".format(root_path))
        self._root_path = root_path
        self._storage_config = storage_config
        self._shear_root_path = root_shear
        self._shear_cat_config = shear_cat_config
        self._cache = cache

    def list_analysis_product_names(self):
        temp = os.listdir(self_root_path)
        return [path for path in temp if os.path.isdir(path)]

    def add_analysis_product_storage(self,name):
        """
        Should probably be using some kind of secondary list where we compare 
        the cosmology collection to allowed cosmology collections.
        """
        if name in self:
            raise ValueError("{:s} already exists".format(name))
        # create the directory
        new_path = os.path.join(self._root_path,name)
        os.mkdir(new_path)
        temp = _load_full_storage_collection(self._root_path,
                                             self._storage_config)
        if self._cache is not None:
            self._cache[name] = temp
        return temp

    def get_analysis_product_storage(self,name):
        if self._cache is not None:
            if name in self._cache:
                return self._cache[name]

        if name not in self:
            raise ValueError("Does not contain {:s}".format(name))
        
        path = os.path.join(self._root_path,name)
        temp = _load_full_storage_collection(self._root_path,
                                             self._storage_config)
        if self._cache is not None:
            self._cache[name] = temp
        return temp

    def get_shear_cat_loader(self,name):
        dir_path = self._shear_root_path.format(name)
        if not os.path.isdir(os.path.normpath(dir_path)):
            raise ValueError(("{:s} does not exist/is not a "
                              "directory").format(dir_path))
        temp = self._shear_cat_config
        return temp.constructCompleteShearCatLoader(dir_path)

    def __contains__(self,name):
        return os.path.isdir(os.path.join(self._root_path,name))


class SampledStorageCollection(CosmologyAnalysisCollection):
    """
    Tracks names of different cosmologies.
    """

class FiducialStorageCollection(CosmologyAnalysisCollection):
    """
    Tracks names of different source configurations.
    """

    def __init__(self,cosmo_name, root, storage_config, root_shear,
                 shear_cat_config, cache = None):
        super(FiducialCosmologyCollection,self).__init__(root, storage_config,
                                                         root_shear,
                                                         shear_cat_config,
                                                         cache = None)
        self._cosmo_name = cosmo_name

    def get_cosmo_name(self):
        return self._cosmo_name

def _get_abs_paths(config,section_option_l,config_file_path):
    """
    Gets absolute paths from a Configuration File.
    """
    
    config_file_dir = os.path.dirname(config_file_path)

    starting_dir = os.getcwd()
    if config_file_dir != '':
        os.chdir(config_file_dir)
    out = []
    for section, option in section_option_l:
        try:
            path = os.path.normpath(config.get(section,option))
        except:
            os.chdir(starting_dir)
            raise
        out.append(os.path.abspath(path))
    if config_file_dir != '':
        os.chdir(starting_dir)
    return out

def _get_paths_cosmo_config_paths(config, config_file_path):
    section_option_l = [("SamplingCosmology","root_dir"),
                        ("SamplingCosmology","storage_config"),
                        ("SamplingCosmology","shear_cat_root_dir"),
                        ("SamplingCosmology","shear_cat_config"),
                        ("FiducialCosmology","root_dir"),
                        ("FiducialCosmology","storage_config"),
                        ("FiducialCosmology","shear_cat_root_dir"),
                        ("FiducialCosmology","shear_cat_config")]
    paths = _get_abs_paths(config,section_option_l,config_file_path)
    sections,options = zip(*section_option_l)
    dict_entries = zip(options,paths)

    sampling_dict = dict(dict_entries[:4])
    fiducial_dict = dict(dict_entries[4:])
    return sampling_dict,fiducial_dict

def _setup_config_parser(config_data):
    #storage_config = ConfigParser.SafeConfigParser()
    #storage_config.read(config_data["storage_config"])
    storage_config = StorageConfig(config_data["storage_config"])

    #shear_cat_config = ConfigParser.SafeConfigParser()
    #shear_cat_config.read(config_data["shear_cat_config"])
    shear_cat_config = ShearCatCollectionLoaderConfig(config_data["shear_cat_config"])
    return storage_config,shear_cat_config

class CosmoCollectionConfigBuilder(object):
    def __init__(self, sampling_cosmology_config,fiducial_cosmology_config):
        self._sampled_cosmology_config = sampling_cosmology_config

        self._fiducial_cosmology_config = fiducial_cosmology_config

    @classmethod
    def from_config_file(cls,fname):
        config = ConfigParser.SafeConfigParser()
        config.read(fname)
        sampling_dict, fiducial_dict = _get_paths_cosmo_config_paths(config,
                                                                     fname)

        # finally let's add the fiducial name to the fiducial dict
        fiducial_dict['fiducial_name'] = config.get("FiducialCosmology",
                                                    "fiducial_name")
        return cls(sampling_dict,fiducial_dict)


    def get_fiducial_storage_collection(self,include_cache = False,
                                        cache_limit=None):
        config_data = self._fiducial_cosmology_config
        storage_config, shear_cat_config = _setup_config_parser(config_data)
        cache = None
        if include_cache:
            raise RuntimeError()

        return FiducialStorageCollection(config_data["fiducial_name"],
                                         config_data['root_dir'],
                                         storage_config,
                                         config_data["shear_cat_root_dir"],
                                         shear_cat_config,
                                         cache)

    def get_sampled_storage_collection(self,include_cache = False,
                                       cache_limit = None):
        config_data = self._sampled_cosmology_config
        storage_config, shear_cat_config = _setup_config_parser(config_data)

        cache = None
        if include_cache:
            raise RuntimeError()

        return SampledStorageCollection(config_data['root_dir'],
                                        storage_config,
                                        config_data["shear_cat_root_dir"],
                                        shear_cat_config,
                                        cache)
