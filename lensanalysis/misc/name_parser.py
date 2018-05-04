from collections import Mapping

from enum_definitions import Descriptor, descriptor_mapping, \
    analysis_object_descriptors, single_descriptor_flags

def _parse_name(name, analysis_names):
    """
    identifies the object type and appropriate descriptors for a name.

    Parameters
    ----------
    name : str
        A name to be parsed. All descriptors in the name are to be separated by 
        an underscore.
    analysis_names : list
        A list of all of the analysis objects.

    Returns
    -------
    descriptors : tuple
        Contains descriptor enums taken from the name.
    object_name : str
        Contains the analysis object name.
    """
    for ref_name in analysis_names:
        ref_length = len(ref_name)
        if ref_name == name[-ref_length:]:
            object_name = ref_name
            break
    else:
        raise ValueError(("{:s} matches no known analysis object "
                          "name").format(name))
    
    if ref_length == len(name):
        return (), name

    len_ref_split = len(ref_name.split('_'))
    name_prefixes = name.split('_')[:-len_ref_split]

    descriptor_l = []
    for prefix in name_prefixes:
        for descriptor_name, descriptor in single_descriptor_flags:
            aliases = descriptor_mapping[descriptor]
            if prefix in aliases:
                # this means that we need to add the descriptor to the
                # descriptor list
                if descriptor in descriptor_l:
                    raise ValueError(("{:s} descriptor specified multiple times"
                                      " in {:s}.").format(descriptor_name,
                                                          name))
                descriptor_l.append(descriptor)
    return tuple(descriptor_l), object_name


class NameParser(object):
    """
    Identifies the object type and appropriate descriptors from a name.

    Parameters
    ----------
    analysis_names : list
        A list of all of the analysis objects.
    """
    def __init__(self,analysis_names = analysis_object_descriptors.keys()):
        self._analysis_names = analysis_names

    def parse_name(self,name):
        return _parse_name(name, self._analysis_names)

def _check_valid_descriptors(descriptors, object_name, mapping,
                             must_contain = True):
    if object_name in mapping or must_contain:
        valid_descriptors = mapping[object_name]
        for elem in descriptors:
            if elem not in valid_descriptors:
                return elem
    return None
            
class SafeNameParser(NameParser):
    """
    Identifies the object type and allowed descriptors from a name.

    The difference between SafeNameParser and NameParser is that name parser 
    checks to see if descriptors are allowed or not

    Parameters
    ----------
    analysis_descr_map : Mapping
        A mapping of all valid analysis objects to sequences of the allowed 
        members of the Descriptors.
    extra_requirements : Mapping, optional
        An optional mapping of a subset of all analysis objects to allowed 
        sequences. These are further restrictions on allowed descriptors to 
        the ones allowed by analysis_descriptor mapping. If this is None, 
        then no optional requirements are checked. Default is None
    """

    def __init__(self,analysis_descr_map = analysis_object_descriptors,
                 extra_requirements = None):
        assert isinstance(analysis_descr_map, Mapping)
        self._analysis_descr_map = analysis_descr_map
        self.set_extra_requirements(extra_requirements)

    def get_extra_requirements(self):
        return self._extra_requirements

    def set_extra_requirements(self,extra_requirements):
        assert (isinstance(extra_requirements,Mapping)
                or extra_requirements is None)
        self._extra_requirements = extra_requirements

    def remove_extra_requirements(self):
        self._extra_requirments = None

    def parse_name(self,name):
        parsed = _parse_name(name, self._analysis_descr_map.keys())

        if out[0] == ():
            with_full_flags = (Descriptor.none,out[1])
        else:
            with_full_flags = parsed

        invalid_descr =_check_valid_descriptors(with_full_flags[0],
                                                with_full_flags[1],
                                                self._analysis_descr_map,
                                                must_contain = True)
        extra_requirements = self.get_extra_requirements()
        if invalid_descr is not None:
            raise ValueError(("{:s} is not a valid descriptor of "
                              "{:s}.").format(str(invalid_descr),
                                              with_full_flags[1]))
        if extra_requirements is not None:
            invalid_descr =_check_valid_descriptors(with_full_flags[0],
                                                    with_full_flags[1],
                                                    extra_requirements,
                                                    must_contain = False)
            if invalid_descr is not None:
                raise ValueError(("{:s} is not a valid descriptor of "
                                  "{:s}.").format(str(invalid_descr),
                                                  with_full_flags[1]))
        return parsed
