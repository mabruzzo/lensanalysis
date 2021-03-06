import ConfigParser
from collections import OrderedDict
import os
import os.path

from .storage_config import StorageConfig, ShearCatCollectionLoaderConfig
from .photoz_config import PseudoPhotozConfig

from ..misc.analysis_collection import UniformAnalysisProductCollection
from ..misc.enum_definitions import Descriptor



def _load_single_storage_collection(path, storage_config_method, descriptions,
                                    tomo_descriptor = False):
    """
    This is probably not the way to do this.
    """
    if tomo_descriptor:
        descriptions = descriptions | Descriptor.tomo


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
                #sanity check
                raise RuntimeError("IT SHOULDN'T BE POSSIBLE FOR THIS TO GET "
                                   "CALLED")
        else:
            raise ValueError("{:s} does not exist/ is not a directory")
    return storage

def _load_conv_storage_collection(path, analysis_collection,
                                  storage_config, save_config=None,
                                  tomo=False):
    """
    Handles the setting up of loading ConvergenceMap Storage Objects.
    """

    kappa_save_config = None
    if tomo:
        ksc = analysis_collection.tomo_conv_map
        if save_config is not None:
            kappa_save_config = save_config.tomo_conv_map
        tomo_prefix = 'Tomographic '
    else:
        ksc = analysis_collection.conv_map
        if save_config is not None:
            kappa_save_config = save_config.conv_map
        tomo_prefix = ''
    no_save_config = (kappa_save_config is None)


    storage_method = storage_config.convergence_map_collection_storage
    if no_save_config or kappa_save_config.noiseless_map:
        ksc.noiseless_map = _load_single_storage_collection(path,
                                                            storage_method,
                                                            Descriptor.none,
                                                            tomo)
        if (not no_save_config) and ksc.noiseless_map is None:
            raise ValueError(("{:s}Noiseless convergence map storage options "
                              "not specified").format(tomo_prefix))

    if no_save_config or kappa_save_config.noisy_map:
        ksc.noisy_map = _load_single_storage_collection(path,
                                                        storage_method,
                                                        Descriptor.noisy,
                                                        tomo)
        if (not no_save_config) and ksc.noisy_map is None:
            raise ValueError(("{:s}Noisy convergence map storage options not "
                              "specified").format(tomo_prefix))

    if no_save_config or kappa_save_config.smoothed_map:
        ksc.smoothed_map = _load_single_storage_collection(path, storage_method,
                                                       Descriptor.smoothed,
                                                       tomo)
        if (not no_save_config) and ksc.smoothed_map is None:
            raise ValueError(("{:s}Smoothed convergence map storage options "
                              "not specified").format(tomo_prefix))

    if no_save_config or kappa_save_config.smoothed_noisy_map:
        temp = _load_single_storage_collection(path, storage_method,
                                               Descriptor.smoothed_noisy,
                                               tomo)
        ksc.smoothed_noisy_map = temp
        if (not no_save_config) and ksc.smoothed_noisy_map is None:
            raise ValueError(("{:s}Smoothed Noisy convergence map storage "
                              "options not specified").format(tomo_prefix))

def _load_shear_storage_collection(path, analysis_collection,
                                   storage_config, save_config=None,
                                   tomo=False):
    """
    Handles the setting up of loading ShearMap Storage Objects.
    """

    shear_save_config = None
    if tomo:
        ssc = analysis_collection.tomo_shear_map
        if save_config is not None:
            shear_save_config = save_config.tomo_shear_map
            tomo_prefix ='Tomographic '
    else:
        ssc = analysis_collection.shear_map
        if save_config is not None:
            shear_save_config = save_config.shear_map
            tomo_prefix = ''
    no_save_config = (shear_save_config is None)

    storage_method = storage_config.shear_map_collection_storage
    if no_save_config or shear_save_config.noiseless_map:
        ssc.noiseless_map = _load_single_storage_collection(path,
                                                            storage_method,
                                                            Descriptor.none,
                                                            tomo)
        if (not no_save_config) and ssc.noiseless_map is None:
            raise ValueError(("{:s}Noiseless shear map storage options "
                              "not specified").format(tomo_prefix))

    if no_save_config or shear_save_config.noisy_map:
        ssc.noisy_map = _load_single_storage_collection(path,
                                                        storage_method,
                                                        Descriptor.noisy,
                                                        tomo)
        if (not no_save_config) and ssc.noisy_map is None:
            raise ValueError(("{:s}Noisy shear map storage options not "
                              "specified").format(tomo_prefix))

def _load_feature_product_storage_collection(path, analysis_collection,
                                             storage_config, save_config=None):
    """
    Handles the setting up of loading FeatureProduct Storage Objects.
    """

    fp_save_config = None
    fpc = analysis_collection.feature_products
    if save_config is not None:
        fp_save_config = save_config.feature_products

    no_save_config = (fp_save_config is None)
    storage_method = storage_config.peak_loc_collection_storage
    if no_save_config or fp_save_config.peak_locations:
    
        fpc.peak_locations = _load_single_storage_collection(path,
                                                             storage_method,
                                                             Descriptor.none)
        if (not no_save_config) and fpc.peak_locations is None:
            raise ValueError(("Peak Location storage options not "
                              "specified"))
    if no_save_config or fp_save_config.tomo_peak_locations:
        temp = _load_single_storage_collection(path, storage_method,
                                               Descriptor.tomo)
        fpc.tomo_peak_locations = temp
        if (not no_save_config) and fpc.tomo_peak_locations is None:
            raise ValueError(("Tomographic Peak Location storage options not "
                              "specified"))

    storage_method = storage_config.peak_counts_collection_storage
    if no_save_config or fp_save_config.peak_counts:
        fpc.peak_counts = _load_single_storage_collection(path,
                                                          storage_method,
                                                          Descriptor.none)
        if (not no_save_config) and fpc.peak_counts is None:
            raise ValueError(("Peak Counts storage options not "
                              "specified"))

    if no_save_config or fp_save_config.tomo_peak_counts:
        fpc.tomo_peak_counts = _load_single_storage_collection(path,
                                                               storage_method,
                                                               Descriptor.tomo)
        if (not no_save_config) and fpc.tomo_peak_counts is None:
            raise ValueError(("Tomographic Peak Counts storage options not "
                              "specified"))

    storage_method = storage_config.power_spectrum_collection_storage
    if no_save_config or fp_save_config.tomo_power_spectrum:
        fpc.tomo_power_spectrum \
            = _load_single_storage_collection(path, storage_method,
                                              Descriptor.tomo)
        if (not no_save_config) and fpc.tomo_power_spectrum is None:
            raise ValueError(("Tomographic Power Spectrum storage options not "
                              "specified"))

    
def _load_full_storage_collection(path,storage_config,save_config = None):
    """
    Create the storage collection as you go.

    For now let's just build everything we can. Probably not a great idea.

    Come back and update after allowing for tomography.
    """
    # the storage collection is uniform in the sense that all
    # AnalysisProductObjects must be None or instances of (virtual) subclasses
    # of CollectionStorage
    analysis_collection = UniformAnalysisProductCollection.factory()

    # kappa storage collection
    ksc = analysis_collection.conv_map
    _load_conv_storage_collection(path, analysis_collection, storage_config,
                                  save_config=save_config, tomo=False)

    # tomographic convergence maps
    _load_conv_storage_collection(path, analysis_collection, storage_config,
                                  save_config=save_config, tomo=True)

    # shear storage collection
    _load_shear_storage_collection(path, analysis_collection, storage_config,
                                   save_config=save_config, tomo=False)
    # tomographic shear storage collection
    _load_shear_storage_collection(path, analysis_collection, storage_config,
                                   save_config=save_config, tomo=True)
    
    
    _load_feature_product_storage_collection(path, analysis_collection,
                                             storage_config,
                                             save_config=save_config)

    
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

    def get_num_tomo_bin(self):
        """
        Get the number of tomographic bins.
        This number is independent of whether or not any tomgraphic data will 
        or can be saved.
        """
        return self._storage_config.num_tomo_bins()

    def list_analysis_product_names(self):
        root_path = self._root_path
        temp = os.listdir(root_path)
        return [path for path in temp if os.path.isdir(os.path.join(root_path,
                                                                    path))]

    def add_analysis_product_storage(self,name,save_config=None):
        """
        Should probably be using some kind of secondary list where we compare 
        the cosmology collection to allowed cosmology collections.
        """
        if name in self:
            raise ValueError("{:s} already exists".format(name))
        # create the directory
        new_path = os.path.join(self._root_path,name)
        os.mkdir(new_path)

        temp = _load_full_storage_collection(new_path,
                                             self._storage_config,
                                             save_config)
        if self._cache is not None:
            self._cache[name] = temp
        return temp

    def get_analysis_product_storage(self,name,save_config=None):
        if self._cache is not None:
            if name in self._cache:
                return self._cache[name]

        if name not in self:
            raise ValueError("Does not contain {:s}".format(name))

        path = os.path.join(self._root_path,name)
        temp = _load_full_storage_collection(path,
                                             self._storage_config,
                                             save_config)
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
                 shear_cat_config, ppz_config, cache = None):
        super(FiducialStorageCollection,self).__init__(root, storage_config,
                                                         root_shear,
                                                         shear_cat_config,
                                                         cache = cache)
        self._cosmo_name = cosmo_name
        self.ppz_config = ppz_config

    def get_cosmo_name(self):
        return self._cosmo_name

    def get_ppz_shear_cat_loader(self,name):
        assert name[-4:] == '_ppz'
        fid_name = self.ppz_config.get_fid_name(name)
        if fid_name is None:
            raise ValueError("No fid_name has been specified for "
                             "{:s}".format(name))
        return self.get_shear_cat_loader(fid_name)
    
def get_abs_paths(config,section_option_l,config_file_path):
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

def _get_cosmo_config_paths(config, config_file_path, section_option_l):
    paths = get_abs_paths(config,section_option_l,config_file_path)
    sections,options = zip(*section_option_l)
    dict_entries = zip(options,paths)
    return dict(dict_entries)

def _get_sampled_cosmo_config_paths(config, config_file_path):
    section_option_l = [("SamplingCosmology","root_dir"),
                        ("SamplingCosmology","storage_config"),
                        ("SamplingCosmology","shear_cat_root_dir"),
                        ("SamplingCosmology","shear_cat_config")]
    return _get_cosmo_config_paths(config, config_file_path, section_option_l)

def _get_fid_cosmo_config_paths(config,config_file_path):
    section_option_l = [("FiducialCosmology","root_dir"),
                        ("FiducialCosmology","storage_config"),
                        ("FiducialCosmology","shear_cat_root_dir"),
                        ("FiducialCosmology","shear_cat_config")]
    photoz = False
    if config.has_option("FiducialCosmology","photoz_config"):
        photoz = True
        section_option_l.append(("FiducialCosmology","photoz_config"))
    temp = _get_cosmo_config_paths(config, config_file_path, section_option_l)
    if not photoz:
        temp["photoz_config"] = None
    return temp

def _setup_config_parser(config_data):
    storage_config = StorageConfig(config_data["storage_config"])
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

        sampling_dict = fiducial_dict = None
        if config.has_section("SamplingCosmology"):
            sampling_dict = _get_sampled_cosmo_config_paths(config, fname)
        if config.has_section("FiducialCosmology"):
            fiducial_dict = _get_fid_cosmo_config_paths(config,fname)
            # finally let's add the fiducial name to the fiducial dict
            fiducial_dict['fiducial_name'] = config.get("FiducialCosmology",
                                                        "fiducial_name")
        return cls(sampling_dict,fiducial_dict)


    def get_fiducial_storage_collection(self,include_cache = False,
                                        cache_limit=None,
                                        save_config = None):
        if self._fiducial_cosmology_config is None:
            raise ValueError("no fiducial cosmology section was specified.")
        config_data = self._fiducial_cosmology_config
        storage_config, shear_cat_config = _setup_config_parser(config_data)
        cache = None
        if include_cache:
            raise RuntimeError()
        if self._fiducial_cosmology_config["photoz_config"] is None:
            ppz_config = None
        else:
            temp = ConfigParser.SafeConfigParser()
            temp.read(self._fiducial_cosmology_config["photoz_config"])
            ppz_config = PseudoPhotozConfig(temp)
        
        return FiducialStorageCollection(config_data["fiducial_name"],
                                         config_data['root_dir'],
                                         storage_config,
                                         config_data["shear_cat_root_dir"],
                                         shear_cat_config,
                                         ppz_config,
                                         cache)

    def get_sampled_storage_collection(self,include_cache = False,
                                       cache_limit = None):
        if self._sampled_cosmology_config is None:
            raise ValueError("no sampled cosmology section was specified.")
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
