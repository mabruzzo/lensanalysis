from abc import ABCMeta, abstractproperty, abstractmethod
from collections import Sequence
import os
import os.path
import warnings

import numpy as np
import astropy.units as u

from lenstools import ShearMap, ConvergenceMap
from lenstools.catalog import ShearCatalog

from fname_formatter import AbstractFnameFormatter
from feature_object import PeakLocations, PeakCounts, TomoPowerSpectra


def write_peak_loc_npy(fname,peak_loc):
    heights = peak_loc.heights
    positions = (peak_loc.locations).to(u.degree).value
    array = np.column_stack((heights,
                             positions[:,0],
                             positions[:,1]))
    np.save(fname,array)

def load_peak_loc_npy(fname):
    array = np.load(fname)
    heights = array[:,0]
    locations = array[:,1:] *u.degree
    return PeakLocations(heights = heights,locations = locations)

def write_peak_count_npy(fname,peak_count):
    counts = peak_count.counts
    np.save(fname,counts)

def load_peak_count_npy(fname):
    return PeakCounts(np.load(fname))

def write_tomo_power_spectra_npy(fname,tomo_power_spectra):
    np.save(fname,tomo_power_spectra)

def read_tomo_power_spectra_npy(fname):
    power = np.load(fname)
    num_spectra = power.shape[0]
    num_tomo_bins = (int(sqrt(8*triangle_num+1))-1)//2

    return TomoPowerSpectra(power,num_tomo_bins)

class CollectionSaver(object):
    __metaclass__ = ABCMeta
    """
    Abstract Base Class for saving an object.
    """
    @abstractmethod
    def save(collection,collection_id):
        pass

    @classmethod
    def __subclasshook__(cls,C):
        if cls is CollectionSaver:
            if any("save"  in B.__dict__ for B in C.__mro__):
                return True
        return NotImplemented

class CollectionLoader(object):
    __metaclass__ = ABCMeta
    """
    Abstract Base Class for saving an object.
    """
    @abstractmethod
    def load(collection_id):
        pass

    @classmethod
    def __subclasshook__(cls,C):
        if cls is CollectionLoader:
            if any("load"  in B.__dict__ for B in C.__mro__):
                return True
        return NotImplemented


class CollectionStorage(object):
    __metaclass__ = ABCMeta
    """
    Abstract Base Class that represents storage of an object.

    Needs to be able to read and save storage.
    """
    @abstractmethod
    def save(collection, collection_id):
        pass

    @abstractmethod
    def load(collection_id):
        pass

    @abstractmethod
    def delete(collection_id):
        pass

    @abstractmethod
    def __contains__(collection_id):
        """
        If only part of the collection is contained, then this returns False.
        """
        pass

    @classmethod
    def __subclasshook__(cls,C):
        if cls is CollectionStorage:
            if any("load"  in B.__dict__ for B in C.__mro__):
                if any("save" in B.__dict__ for B in C.__mro__):
                    if any("delete" in B.__dict__ for B in C.__mro__):
                        if any("__contains__" in B.__dict__ for B in C.__mro__):
                            return True
        return NotImplemented

def _check_field_mapping(fname_formatter, num_elements, field_mapping):
    map_keys = field_mapping.keys()
    fields = fname_formatter.fields
    if len(map_keys) not in [1,2]:
        raise ValueError(("map_keys must contain {:d} "
                          "entries.").format(len(fields)))
    if "collection_id" not in map_keys:
        raise ValueError("collection_id must be 1 of the keys in "
                         "field_mapping.")

    cid_field = field_mapping["collection_id"]
    if cid_field in fields:
        remaining = [field for field in fields if field != cid_field ]
    else:
        raise ValueError("collection_id does not map to a valid field.")

    if len(remaining) > 1:
        raise ValueError("fname_formatter has too many fields.")


    if num_elements == 1:
        fields = fname_formatter.fields
        if len(remaining) == 0:
            if len(map_keys) == 1:
                return None

    # we can still have 2 fields if there is only 1 element
    if len(remaining) != 1:
        raise ValueError("fname_formatter should have 2 fields.")
    if "element_id" not in map_keys:
        raise ValueError("element_id must be a key in field_mapping.")
    if field_mapping["element_id"] != remaining[0]:
        if field_mapping["element_id"] == field_mapping["collection_id"]:
            raise ValueError("element_id and collection_id must be mapped to "
                             "different fields.")
        else:
            raise ValueError(("element_id must be mapped to "
                              "{:s}.").format(remaining[0]))


class _BaseFileGroupCollection(object):
    """
    Base class for representing elements of collections as individual files

    Parameters
    ----------
    fname_formatter : AbstractFnameFormatter
        An object with the method format_fname that formats the file name
    root_dir : str
        The root directory in which all of the files are stored
    num_elements : int
        The number of entries of all collections stored here
    field_mapping : dict
        A dictionary that maps "collection_id" and "element_id" to the fields 
        of property of fname_formatter. If num_elements is 1, and there is only 
        1 field in fname_formatter, then this does not need to contain an 
        entry for "element_id."
    """
    def __init__(self, fname_formatter, root_dir, num_elements, field_mapping):
        if isinstance(fname_formatter,AbstractFnameFormatter):
            self._fname_formatter = fname_formatter
        else:
            raise TypeError("fname_formatter must be a (virtual) subclass of "
                            "the abc, AbstractFnameFormatter")

        self.root_dir = root_dir
        if not isinstance(num_elements,int):
            raise TypeError("num_elements must be an int")
        elif num_elements<=0:
            raise ValueError("num_elements must be at leat 1")

        _check_field_mapping(fname_formatter, num_elements, field_mapping)

        self._num_elements = num_elements
        self._cid_field = field_mapping["collection_id"]
        if "element_id" in field_mapping:
            self._eid_field = field_mapping["element_id"]
        else:
            self._eid_field = None

    def _format_fname(self,cid,eid):
        """
        Private helper method that finds the fname for a given collection 
        element.
        """
        if self._num_elements == 1 and self._eid_field is None:
            fname = self._fname_formatter.format_fname(**{self._cid_field:cid})
        else:
            fname = self._fname_formatter.format_fname(**{self._cid_field:cid,
                                                          self._eid_field:eid})
        return '/'.join((self.root_dir,fname))

def _construct_subdir(root_dir,subdir_l):
    """
    Helper function to ensure that all necessary subdirectories exist.
    """
    for subdir,subsequent_subdir in subdir_l:
        subdir_path = os.path.join(root_dir,subdir)
        if not os.path.isdir(subdir_path):
            if os.path.exists(subdir_path):
                raise RuntimeError(("{} already exists but is not a "
                                    "directory.").format(subdir_path))
            os.mkdir(subdir_path)
        if subsequent_subdir is not None:
            _construct_subdir(subdir_path,subsequent_subdir)

class FileGroupCollectionStorage(_BaseFileGroupCollection):
    """
    Base class for representing elements of collections as individual files.

    It uses 2 class variables to contain callables that read and write 
    individual elements in the collections that this is storage for. These are 
    to be overwritten in subclasses.

    Parameters
    ----------
    fname_formatter : AbstractFnameFormatter
        An object with the method format_fname that formats the file name
    root_dir : str
        The root directory in which all of the files are stored
    num_elements : int
        The number of entries of all collections stored here
    field_mapping : dict
        A dictionary that maps "collection_id" and "element_id" to the fields 
        of property of fname_formatter. If num_elements is 1, and there is only 
        1 field in fname_formatter, then this does not need to contain an 
        entry for "element_id."
    """

    element_writer = None
    element_loader = None

    @property
    def num_elements(self):
        return self._num_elements

    @property
    def fname_formatter(self):
        return self._fname_formatter

    def save(self,collection, collection_id):
        if self._num_elements != len(collection):
            raise ValueError(("collection must have {:d} "
                              "elements.").format(self._num_elements))
        for i,elem in enumerate(collection):
            fname = self._format_fname(collection_id,i+1)
            self.element_writer(fname,elem)

    def load(self,collection_id):
        if collection_id not in self:
            raise ValueError(("collection_id {:d} is not contained within" 
                              "storage").format(collection_id))

        out = []
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            out.append(self.element_loader(fname))
        return out

    def delete(self,collection_id):
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            if os.path.exists(fname):
                os.remove(fname)

    def __contains__(self,collection_id):
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            if os.path.exists(fname):
                if os.path.isfile(fname):
                    continue
            return False
        return True

    def construct_subdirectories(self, start_realization,stop_realization):
        """
        Ensures that all subdirectories exist.

        While we would probably be better off tring to create the subdir as 
        they are needed (instead of all at once like this method does), this 
        would require an overhaul of the classes so that when the program is 
        run in parallel, there is no possibility of conflicts.

        Parameters
        ----------
        start_realization : int
            The minimum realization that may need to be saved
        stop_realization : int
            stop_realization-1 is the last realization that may need to be 
            saved.
        """

        subdir_l = self.fname_formatter.determine_subdir(start_realization,
                                                         stop_realization-1)
        if subdir_l is None:
            return
        _construct_subdir(self.root_dir,subdir_l)

_shear_map_writer = lambda fname, shear_map : shear_map.save(fname)
class ShearMapCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for shear map collections
    """
    # note that element_writer can't be linked directly to the lambda function
    # unless it takes an additional argument, self, in which case it will be 
    # an instance method
    element_writer = staticmethod(_shear_map_writer)
    element_loader = staticmethod(ShearMap.load)

_conv_map_writer = lambda fname, conv_map : conv_map.save(fname)
class ConvergenceMapCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for convergence map collections
    """

    element_writer = staticmethod(_conv_map_writer)
    element_loader = staticmethod(ConvergenceMap.load)


class PeakLocCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for peak location collections
    """

    element_writer = staticmethod(write_peak_loc_npy)
    element_loader = staticmethod(load_peak_loc_npy)


class PeakCountCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for peak count collections
    """

    element_writer = staticmethod(write_peak_count_npy)
    element_loader = staticmethod(load_peak_count_npy)


class FullShearCatFGLoader(_BaseFileGroupCollection):
    """
    Loads a full shear catalog from files representing individual catalogs.

    Parameters
    ----------
    fname_formatter : AbstractFnameFormatter
        An object with the method format_fname that formats the file name
    root_dir : str
        The root directory in which all of the files are stored
    num_elements : int
        The number of entries of all collections stored here
    field_mapping : dict
        A dictionary that maps "collection_id" and "element_id" to the fields 
        of property of fname_formatter. If num_elements is 1, and there is only 
        1 field in fname_formatter, then this does not need to contain an 
        entry for "element_id."
    position_fname_formatter : AbstractFnameFormatter
        An object with the method format_fname that formats the file name of 
        the file that track position information of a file.
    """
    def __init__(self, fname_formatter, root_dir, num_elements,
                 field_mapping, position_fname_formatter):
        super(FullShearCatFGLoader,self).__init__(fname_formatter, root_dir,
                                                  num_elements, field_mapping)

        if isinstance(position_fname_formatter,AbstractFnameFormatter):
            self._pos_fname_formatter = position_fname_formatter
        else:
            raise TypeError("position_fname_formatter must be a (virtual) "
                            "subclass of the abc, AbstractFnameFormatter")
        # the whole idea is that the position_fname_formatter is independent of
        # realization - it should only depend on relization.
        fields = position_fname_formatter.fields
        if self._eid_field is None:
            if len(fields) != 0:
                raise ValueError("position_fname_formatter should have no "
                                 "fields")
        if len(fields) >1:
            raise ValueError("position_fname_formatter should have 1 field.")
        if fields[0] != self._eid_field:
            raise ValueError(("position_fname_formatter should have 1 field: "
                              "{:s}").format(self._eid_field))

    def _format_pos_fname(self,eid):
        """
        Private helper method that finds the fname for a given collection 
        element.
        """
        if self._num_elements == 1 and self._eid_field is None:
            fname = self._pos_fname_formatter.format_fname()
        else:
            kwargs = {self._eid_field : eid}
            fname = self._pos_fname_formatter.format_fname(**kwargs)
        return '/'.join((self.root_dir,fname))
        
    def load(self,collection_id):
        out = []
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            pos_fname = self._format_pos_fname(i+1)
            try:
                temp = ShearCatalog.readall([fname],[pos_fname])
            except:
                print "Could not load Shear catalog"
                raise
            out.append(temp)
        return out

class SingleFileCollectionStorage(FileGroupCollectionStorage):
    """
    Base class for representing elements of collections in an individual file.

    It uses 2 class variables to contain callables that read and write 
    individual elements in the collections that this is storage for. These are 
    to be overwritten in subclasses.

    Parameters
    ----------
    fname_formatter : AbstractFnameFormatter
        An object with the method format_fname that formats the file name
    root_dir : str
        The root directory in which all of the files are stored
    num_elements : int
        The number of entries of all collections stored here
    field_mapping : dict
        A dictionary that maps "collection_id" to the field 
        of property of fname_formatter.

    Notes
    ------
    Although, we have chosen to make this class a subclass of 
    FileGroupCollectionStorage for now, it should not be. We get away with this
    by internally using methods that assume that assume that we are always 
    looking for a collection of 1 file. The class hierarchy should 
    be revisited and refactored
    """

    def __init__(self, fname_formatter, root_dir, field_mapping):
        if isinstance(fname_formatter,AbstractFnameFormatter):
            self._fname_formatter = fname_formatter
        else:
            raise TypeError("fname_formatter must be a (virtual) subclass of "
                            "the abc, AbstractFnameFormatter")

        self.root_dir = root_dir
        if not isinstance(num_elements,int):
            raise TypeError("num_elements must be an int")
        elif num_elements<=0:
            raise ValueError("num_elements must be at leat 1")

        if "element_id" in field_mapping:
            raise ValueError("SingleFileCollectionStorage instances should not "
                             "have an 'element_id' key in field_mapping")
        _check_field_mapping(fname_formatter, 1, field_mapping)

        self._num_elements = 1
        self._eid_field = None
        self._cid_field = field_mapping["collection_id"]

class TomoPowerSpectraCollectionSFStorage(FileGroupCollectionStorage):
    """
    Single File storage subclass for tomographic power spectra.
    """
    element_writer = staticmethod(write_tomo_power_spectra_npy)
    element_reader = staticmethod(read_tomo_power_spectra_npy)
