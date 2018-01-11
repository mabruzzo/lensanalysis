from abc import ABCMeta, abstractproperty, abstractmethod
from collections import Sequence
import os
import os.path
import warnings

import numpy as np
import astropy.units as u

from lenstools import ShearMap, ConvergenceMap
from lenstools.catalog import ShearCatalog

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

    @abstractproperty
    def fields(self):
        """
        list of unique fields.
        """
        pass

    @abstractmethod
    def subdirectories(self,**kwargs):
        """
        Returns the subdirectories within which a file is stored in reverse 
        depth order
        """
        pass

    @abstractmethod
    def format_fname(self, **kwargs):
        """
        Apply the procedural step to the packet.
        """
        pass

    @classmethod
    def __subclasshook__(cls,C):
        """
        Need to update.
        """
        if cls is AbstractFnameFormatter:
            if any("format_fname"  in B.__dict__ for B in C.__mro__):
                return True
        return NotImplemented

def _check_template_fields(fname_template,fields):
    if not isinstance(fname_template,basestring):
        raise ValueError("fname_template must be a str")

    if isinstance(fields,basestring):
        fields = (fields,)
    elif isinstance(fields,Sequence):
        l = []
        for elem in fields:
            if not isinstance(elem,basestring):
                raise ValueError("fields must either be an str or a "
                                 "sequence of str")
            l.append(elem)
        fields = tuple(l)
    else:
        raise ValueError("fields must either be an str or a "
                         "sequence of str")

    warnings.warn("Need to implement check that an adequate number of fields "
                  "have been provided.", RuntimeWarning)
    
    return fields

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
        self._fields = _check_template_fields(fname_template,fields)
        self._fname_template = fname_template
        self._unique_fields = list(set(self._fields))

    @property
    def fname_template(self):
        return self._fname_template

    @property
    def fields(self):
        return self._unique_fields

    def subdirectories(self,**kwargs):
        return []

    def format_fname(self,**kwargs):
        """
        Determines a filename given arguments.

        Either specify all arguments with positional arguments (in the order of 
        the listed fields) or give all arguments as kwargs. Currently, we do 
        not accept a mix.
        """

        l = []
        for field in self._fields:
            if field in kwargs:
                l.append(kwargs[field])
            else:
                raise ValueError("Didn't specify a value for %s" % field)
        args = tuple(l)

        return self._fname_template.format(*args)

    
class FnameFormatterDecorator(AbstractFnameFormatter):
    __metaclass__ = ABCMeta
    
    @abstractproperty
    def wrapped_formatter(self,value):
        pass

class SubdirectoryFnameFormatter(FnameFormatterDecorator):

    __metaclass__ = ABCMeta

    @property
    def wrapped_formatter(self):
        return self._wrapped_formatter

    @property
    def fields(self):
        return self._unique_fields

    @abstractmethod
    def _determine_dir_name(self,**kwargs):
        pass

    def subdirectories(self,**kwargs):
        temp = self._wrapped_formatter.subdirectories(**kwargs)
        temp.append(self._determine_dir_name(**kwargs))
        return kwargs

    def format_fname(self,**kwargs):
        temp = os.path.join(self._determine_dir_name(**kwargs),
                            self.wrapped_formatter.format_fname(**kwargs))
        return os.path.normpath(temp)

class RealizationBinnedSubdirectoryFormatter(SubdirectoryFnameFormatter):
    """
    Writes the paths for subdirectories that bin realizations.

    Examples include 1-288, 289-576, etc.

    The realization numbers are assumed to be one-indexed.

    Parameters
    ----------
    directory_template :str
        New style string formatted string for template.
    template_fields : sequence or str
        Names of inserted fields for template. The only allowed fields are 
        'min' (minimum included realization number), 'max' (maximum included 
        realization number), and 'mid' (middle number of included realizations).
    realization_per_dir : int
        Nominal number of realizations per directory.
    wrapped : AbstractFnameFormatter
        Instance of a (virtual) subclass of AbstractFnameFormatter.
    realization_field : str
        The field name associated with realization
    """
     
    def __init__(self, directory_template, template_fields, realization_per_dir,
                 wrapped, realization_field):
        template_fields = _check_template_fields(directory_template,
                                                 template_fields)
        for elem in template_fields:
            if elem not in ["max", "min", "mid"]:
                raise ValueError(("{:s} is an invalid realization binned "
                                  "subdirectory template field.").format(elem))
        self._template_fields = template_fields

        self._directory_template = directory_template

        if isinstance(realization_per_dir,int):
            if realization_per_dir>0:
                self._realization_per_dir = realization_per_dir
            else:
                raise ValueError("realization_per_dir must be at least 1.")
        else:
            raise TypeError("realization_per_dir must be an int")

        if isinstance(realization_field,str):
            self._realization_field = realization_field
        else:
            raise ValueError("realization_field must be a str")

        self._unique_fields = list(set(wrapped.fields +
                                       [self._realization_field]))

        self._wrapped_formatter = wrapped

    def _determine_dir_name(self,**kwargs):
        realization = kwargs[self._realization_field]
        assert isinstance(realization,int)

        temp = realization-1
        bin_num = (realization-1)/self._realization_per_dir
        bin_min = (self._realization_per_dir * bin_num)+1
        bin_max = bin_min + self._realization_per_dir -1

        args = []
        for elem in self._template_fields:
            if elem == 'min':
                args.append(bin_min)
            elif elem == 'max':
                args.append(bin_max)
            elif elem == 'mid':
                args.append((bin_min+bin_max)/2)
        temp = self._directory_template.format(*args)
        return temp

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
            # except Writing error, we can try to create the necessary
            # directories
            # can get the directories needed from : self._fname_formatter
            raise RuntimeError()

    def load(self,collection_id):
        if collection_id not in self:
            raise ValueError("collection_id is not contained within storage")

        out = []
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            out.append(self.element_load(fname))
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
        
    def loader(self,collection_id):
        out = []
        for i in range(self._num_elements):
            fname = self._format_fname(collection_id,i+1)
            pos_fname = self._format_fname(collection_id,i+1)
            try:
                temp = ShearCatalog.readall([fname],[pos_fname])
            except:
                print "Could not load Shear catalog"
                raise
            out.append(temp)
        return out
        
if __name__ == "__main__":
    fname_formatter = BaseFnameFormatter("WLshear_cat_bin{:d}_{:04d}r.fits",
                                         ("bin_num","realization"))
    print
    print fname_formatter.format_fname(**{"realization":27, "bin_num" : 4})
    print fname_formatter.fields

    full_formatter = RealizationBinnedSubdirectoryFormatter('{:d}-{:d}',
                                                            ('min','max'),
                                                            288,
                                                            fname_formatter,
                                                            'realization')
    print full_formatter.fields
    print full_formatter.format_fname(**{"realization":27, "bin_num" : 4})
    print full_formatter.format_fname(**{"realization":287, "bin_num" : 1})
    print full_formatter.format_fname(**{"realization":288, "bin_num" : 2})
    print full_formatter.format_fname(**{"realization":289, "bin_num" : 3})
    print full_formatter.format_fname(**{"realization":576, "bin_num" : 5})
    print full_formatter.format_fname(**{"realization":577, "bin_num" : 4})
