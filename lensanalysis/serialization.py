from abc import ABCMeta
from collections import Sequence
import os
import os.path

import numpy as np
from astropy.units as u

from lenstools import ShearMap, ConvergenceMap
from lenstools.shear import ShearCatalog

from feature_object import PeakLocations, PeakCounts


def write_peak_loc_npy(fname,peak_loc):
    heights = peak_loc.heights
    positions = (peak_loc.locations).to(u.degree).value
    array = np.column_stack((heights,
                             positions[:,0].value,
                             positions[:,1].value))
    np.save(fname,array)

def load_peak_loc_npy(fname):
    array = np.load(fname)
    heights = array[:,0]
    locations = array[:,1:]
    return PeakLocations(heights = heights,locations = locations)
    
def write_peak_count_npy(fname,peak_count):
    counts = peak_count.counts
    np.save(fname,counts)

def load_peak_count_npy(fname):
    return PeakCounts(np.load(fname))

class AbstractFnameFormatter(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def format_fname(self, data_object, packet):
        """
        Apply the procedural step to the packet.
        """
        pass

    @classmethod
    def __subclasshook__(cls,C):
        if cls is AbstractFnameFormatter:
            if any("format_fname"  in B.__dict__ for B in C.__mro__):
                return True
        return NotImplemented


class BaseFnameFormatter(AbstractFnameFormatter):
    """
    Object responsible for formatting the file name.

    Note this class always formats the file name component of a path. It 
    usually needs to be joined with the root directory of where the file is 
    stored

    When passing in the fname_template use new string formatting.

    Need to add code to check to see that this is set up correctly.

    Parameters
    ----------
    fname_template : str
        The formatting template of the file. Use new style string formatting.
    fields : sequence or str
        The names of each field being inserted into the template
    """

    def __init__(self,fname_template,fields):
        if isinstance(fname_template,basestring):
            self._fname_template = fname_template
        else:
            raise ValueError("fname_template must be a str")

        if isinstance(fields,basestring):
            self._fields = (fields,)
        elif isinstance(fields,Sequence):
            l = []
            for elem in fields:
                if not isinstance(fields,basestring):
                    raise ValueError("fields must either be an str or a "
                                     "sequence of str")
                l.append(elem)
            self._fields = tuple(l)
        else:
            raise ValueError("fields must either be an str or a "
                             "sequence of str")

    @property
    def fname_template(self):
        return self._fname_template

    @property
    def fields(self):
        return self._fields

    def format_fname(self,**kwargs):
        """
        Determines a filename given arguments.

        Either specify all arguments with positional arguments (in the order of 
        the listed fields) or give all arguments as kwargs. Currently, we do 
        not accept a mix.
        """

        l = []
        for field in self.fields:
            if field in kwargs:
                l.append(kwarg[field])
            else:
                raise ValueError("Didn't specify a value for %s" % field)
        args = tuple(l)

        return self._fname_template.format(*args)


class CollectionStorage(object):
    __metaclass__ = ABCMeta
    """
    Class that represents storage of an object.

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
        remaining = [field if field != cid_field for field in fields]
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
    if len(remaining != 1):
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

    
class FileGroupCollectionStorage(AbstractCollectionStorage):
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

    def __init__(self, fname_formatter, root_dir, num_elements, field_mapping):
        if isinstance(fname_formatter,AbstractFnameFormatter):
            self._fname_formatter = file_formatter
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

    @property
    def num_elements(self):
        return self._num_elements

    @property
    def fname_formatter(self):
        return self._fname_formatter

    def _format_fname(self,cid,eid):
        """
        Private helper method that finds the fname for a given collection 
        element.
        """
        if self._num_elements == 1 and self._eid_field is None:
            fname = self._fname_formatter.format_fname(self._cid_field=cid)
        else:
            fname = self._fname_formatter.format_fname(self._cid_field=cid,
                                                       self._eid_field=eid)
        return '/'.join((self.root_dir,fname))

    def save(collection, collection_id):
        if self._num_elements != len(collection):
            raise ValueError(("collection must have {:d} "
                              "elements.").format(self._num_elements))
        for i,elem in enumerate(collection):
            fname = self._format_fname(collection_id,i+1)
            self.element_writer(fname,elem)

    def load(collection_id):
        if collection_id not in self:
            raise ValueError("collection_id is not contained within storage")

        out = []
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            out.append(self.element_load(fname))
        return out

    def delete(collection_id):
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

class ShearCatCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for shear catalog collections
    """

    element_writer = lambda fname, shear_cat : shear_cat.write(fname)
    element_loader = ShearCatalog.readall

    def __init__(self, fname_formatter, root_dir, num_elements, field_mapping):
        raise NotImplementedError("Need to revisit this because we can "
                                  "constuct shear maps with/with out z info")

class ShearMapCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for shear map collections
    """

    element_writer = lambda fname, shear_map : shear_map.write(fname)
    element_loader = ShearMap.load

class ConvergenceMapCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for convergence map collections
    """

    element_writer = lambda fname, conv_map : conv_map.save(fname)
    element_loader = ConvergenceMap.load

class PeakLocCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for peak location collections
    """

    element_writer = write_peak_loc_npy
    element_loader = load_peak_loc_npy

class PeakCountCollectionFGStorage(FileGroupCollectionStorage):
    """
    File group storage subclass for peak count collections
    """

    element_writer = write_peak_count_npy
    element_loader = load_peak_count_npy
