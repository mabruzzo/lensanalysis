from collections import MutableMapping

from enum_definitions import Descriptor
from name_parser import SafeNameParser

def _check_sequence_contents(sequence, type):
    """
    makes sure that a sequence only contains objects of a given type
    """
    for elem in sequence:
        if not isinstance(elem,type):
            return elem


def _check_name_format(*args):
    """
    Can be a single string, 
    Can be a single tuple: ((descriptors,),object_name)
                           (combined_descriptor, object_name)
    Can be 2 separate arguments: (descriptors,) object_name
                                  combined_descriptor, object_name
    """
    if len(args) == 1:
        arg = args[0]
        if isinstance(arg,basestring):
            return arg
        elif not isinstance(arg,(tuple,list)):
            raise TypeError("valid single arguments can only be a string, "
                            "tuple or list")
        if len(arg) != 2:
            raise ValueError("tuple/list single arguments must have a length "
                             "of 2")
        first,second = arg
    elif len(args)!= 2:
        raise ValueError("Only 1 or 2 arguments are accepted")
    else:
        first,second = args

    if isinstance(first, (tuple,list)):
        result = _check_sequence_contents(first, Descriptor)
        if result is not None:
            raise ValueError(("descriptors specified by a tuple/list, must "
                              "only contain Descriptor flag objects. {:r} "
                              "is not a Descriptor.").format(result))
    elif not isinstance(first, Descriptor):
        raise ValueError("descriptors must be a specified by a Descriptor Flag "
                         "object or a tuple/list, must of Descriptor flag "
                         "objects")

    if not isinstance(second,basestring):
        raise TypeError("object names must be strings")
    return first,second


class NameFormatAliasConversion(object):
    """
    Maps between different name formats.
    """

    def __init__(self,name_parser):
        self.name_parser = name_parser

    def convert_name(self,*args):
        arg = _check_name_format(*args)

        if isinstance(arg,basestring):
            if self.name_parser is None:
                return ValueError("Can't convert strings because no name "
                                  "parser has been specified")
            else:
                return self.convert_name(self.name_parser.parse_name(arg))

        else:
            if isinstance(arg[0], Descriptor):
                descr = arg[0]
            else:
                descr = Descriptor.none
                for descriptor in arg[0]:
                    descr = descr | descriptor
            return descr,arg[1]

class AnalysisObjectMapping(MutableMapping):
    """
    This maps analysis object descriptors to objects.

    Note that the design of having this wrap NameFormatAliasConversion which 
    wraps SafeNameParser doesn't seem to be the best design. It is currently 
    designed in this way because of the way in which I improved on the overall 
    design. We still may want to improve a little more on where these objects 
    get initialized from.
    """
    _default_converter = NameFormatAliasConversion(SafeNameParser())

    def __init__(self,name_format_converter=None):
        if name_format_converter is not None:
            self._name_format_converter = name_format_converter
        else:
            self._name_format_converter =self._default_converter
        self._dict = {}
        
    def __getitem__(self,key):
        temp = self._name_format_converter.convert_name(key)
        return self._dict[temp]

    def __setitem__(self, key, value):
        temp = self._name_format_converter.convert_name(key)
        self._dict[temp] = value

    def __delitem__(self, key):
        temp = self._name_format_converter.convert_name(key)
        self._dict.__delitem__(temp)

    def keys():
        return self._dict.keys()

    def __iter__(self):
        return self._dict.__iter__()

    def __len__(self):
        return self._dict.__len__()

    def iteritems(self):
        return self._dict.iteritems()

    def __repr__(self):
        return "AnalysisObjectMapping.from_mapping({0!r})".format(self._dict)
    
    @classmethod
    def from_sequence(cls, sequence, name_format_converter=None):
        return cls.from_mapping(dict(sequence),name_format_converter)

    @classmethod
    def from_mapping(cls, mapping, name_format_converter=None):
        out = cls(name_format_converter)
        for key,value in mapping.iteritems:
            out[key] = value
        return out

if __name__ == '__main__':
    t = AnalysisObjectMapping()
    t["smooth_conv_map"] = None
    print t["smooth_conv_map"]
    t["shear_map"] = True
    print repr(t)
    
