from collections import Mapping

from analysis_object_mapping import AnalysisObjectMapping
from enum_definitions import Descriptor


class ConvergenceMapProductCollection(object):
    mapping = [(Descriptor.none, "noiseless_map"),
               (Descriptor.noisy, "noisy_map"),
               (Descriptor.smoothed, "smoothed_map"),
               (Descriptor.smoothed_noisy, "smoothed_noisy_map")]

    def __init__(self):
        self.noiseless_map = None
        self.noisy_map = None
        self.smoothed_map = None
        self.smoothed_noisy_map = None

class ShearMapProductCollection(object):
    mapping = [(Descriptor.none, "noiseless_map"),
               (Descriptor.noisy, "noisy_map")]

    def __init__(self):
        self.noiseless_map = None
        self.noisy_map = None

class FeatureProductCollection(object):
    """
    If we add more features, this will be extended.
    """

    mapping = [(Descriptor.none, "peak_loc", "peak_locations"),
               (Descriptor.tomo, "peak_loc", "tomo_peak_locations"),
               (Descriptor.none, "peak_counts", "peak_counts"),
               (Descriptor.tomo, "peak_counts", "tomo_peak_counts")]
    
    def __init__(self):
        self.peak_locations = None
        self.peak_counts = None
        self.tomo_peak_locations = None
        self.tomo_peak_counts = None

class AnalysisProductCollection(object):
    """
    Keeps track of the analysis products.
    """
    conv_map_sub_collection = ConvergenceMapProductCollection
    shear_map_sub_collection = ShearMapProductCollection
    feature_product_sub_collection = FeatureProductCollection

    def __init__(self):
        self.conv_map = None
        self.tomo_conv_map = None
        self.shear_map = None
        self.tomo_shear_map = None
        self.feature_products = None

    @classmethod
    def factory(cls,conv_map_sub_collection = None,
                shear_map_sub_collection = None,
                feature_product_sub_collection = None):
        if conv_map_sub_collection is None:
            conv_map_sub_collection = cls.conv_map_sub_collection
        if shear_map_sub_collection is None:
            shear_map_sub_collection = cls.shear_map_sub_collection
        if feature_product_sub_collection is None:
            feature_product_sub_collection = cls.feature_product_sub_collection

        out = cls()
        collection_class = conv_map_sub_collection
        out.conv_map = collection_class()
        out.tomo_conv_map = collection_class()

        collection_class = shear_map_sub_collection
        out.shear_map = collection_class()
        out.tomo_shear_map = collection_class()

        collection_class = feature_product_sub_collection
        print collection_class
        out.feature_products = collection_class()
        return out


def _set_prop(self,attr_name,value,object_name=None,
              add_tomo_descr=False):
    assert (not (add_tomo_descr and object_name is None))
    
    if getattr(self,attr_name) is None:
        # Update the value of the attribute
        setattr(self,attr_name,value)

        # update the internal map to keep track of how to dynamically access
        # all members
        mapping = getattr(self,"_mapping")

        analysis_object_name = object_name

        for elem in value.mapping:
            if object_name is None:
                assert len(elem) == 3
                descriptor,analysis_object_name,sub_attr_name = elem
            else:
                assert len(elem) == 2
                descriptor,sub_attr_name = elem
                if add_tomo_descr:
                    descriptor = descriptor | Descriptor.tomo

            mapping[(descriptor,analysis_object_name)] = (attr_name,
                                                          sub_attr_name)
    else:
        raise ValueError("Cannot modify the property because it has already "
                         "been initizalized")

def _get_prop(self,attr_name):
    return getattr(self,attr_name)

class UniformAnalysisProductCollection(AnalysisProductCollection,Mapping):
    """
    A subclass of AnalysisProductCollection designed to hold uniform data 
    consistent between all component data structure containers. It has the 
    Mapping interface to allow for dynamic access of component data structure 
    details.

    This is accomplished by making each of the component data structures 
    read_only.

    We dynamically set up properties to minimize code duplication
    """

    def __init__(self):
        self._mapping = AnalysisObjectMapping()
        self._conv_map = None
        self._tomo_conv_map = None
        self._shear_map = None
        self._tomo_shear_map = None
        self._feature_products = None

    conv_map = property(lambda self: _get_prop(self,'_conv_map'),
                        lambda self,value : _set_prop(self,'_conv_map',
                                                      value, 'conv_map'))
    tomo_conv_map = property(lambda self: _get_prop(self,'_tomo_conv_map'),
                             lambda self,value : _set_prop(self,
                                                           '_tomo_conv_map',
                                                           value, 'conv_map',
                                                           True))
    shear_map = property(lambda self: _get_prop(self,'_shear_map'),
                         lambda self,value : _set_prop(self,'_shear_map',
                                                       value, 'shear_map'))
    tomo_shear_map = property(lambda self: _get_prop(self,'_tomo_shear_map'),
                              lambda self,value : _set_prop(self,
                                                            '_tomo_shear_map',
                                                            value, 'shear_map',
                                                            True))
    feature_products=property(lambda self: _get_prop(self,'_feature_products'),
                              lambda self,value:_set_prop(self,
                                                          '_feature_products',
                                                          value))
    
    def __getitem__(self,key):
        first,second = self._mapping[key]
        return getattr(getattr(self,first),second)

    def __setitem__(self,key,value):
        first,second = self._mapping[key]
        setattr(getattr(self,first),second,value)

    def __len__(self):
        return len(self._mapping)

    def __iter__(self):
        return self._mapping.__iter__()

def default_value_UAPC(value):
    """
    Initialize UniformAnalysisProductCollection with default values.
    """
    out = UniformAnalysisProductCollection.factory()
    for key in out:
        out[key] = value
    return out
    
def get_analysis_col_value(analysis_col, descriptors, object_name):
    first,second = _get_attribute_tuple(descriptors, object_name)
    temp = getattr(analysis_col,first)
    return getattr(temp, second)

def _get_attribute_tuple(descriptors, object_name):
    if object_name == "conv_map":
        if Descriptor.tomo in descriptors:
            first = 'tomo_conv_map'
        else:
            first = 'conv_map'

        if Descriptor.noisy in descriptors:
            if Descriptor.smoothed in descriptors:
                second = 'smoothed_noisy_map'
            else:
                second = 'noisy_map'
        elif Descriptor.smoothed in descriptors:
            second = 'smoothed_map'
        else:
            second = 'noiseless_map'

    elif object_name == "shear_map":
        if Descriptor.tomo in descriptors:
            first = 'tomo_shear_map'
        else:
            first = 'shear_map'

        if Descriptor.noisy in descriptors:
            second = 'noisy_map'
        else:
            second = 'noiseless_map'

    elif object_name in ["peak_loc","peak_counts"]:
        first = 'feature_products'

        if object_name == "peak_loc":
            if Descriptor.tomo in descriptors:
                second = "tomo_peak_locations"
            else:
                second = "peak_locations"
        elif Descriptor.tomo in descriptors:
            second = "tomo_peak_counts"
        else:
            second = "peak_counts"
    else:
        raise ValueError()
    return first,second



def set_analysis_col_value(analysis_col, descriptors, object_name, value):
    first,second = _get_attribute_tuple(descriptors, object_name)
    temp = getattr(analysis_col,first)
    setattr(temp, second,value)


