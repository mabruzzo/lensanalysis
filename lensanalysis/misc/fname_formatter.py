from abc import ABCMeta, abstractmethod, abstractproperty
import os.path
import warnings

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
