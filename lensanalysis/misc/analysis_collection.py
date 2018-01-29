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
    def factory(cls):
        out = cls()
        collection_class = cls.conv_map_sub_collection
        out.conv_map = collection_class()
        out.tomo_conv_map = collection_class()

        collection_class = cls.shear_map_sub_collection
        out.shear_map = collection_class()
        out.tomo_shear_map = collection_class()

        collection_class = cls.feature_product_sub_collection
        out.feature_products = collection_class()
        return out
        
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

def get_analysis_col_value(analysis_col, descriptors, object_name):
    first,second = _get_attribute_tuple(descriptors, object_name)
    temp = getattr(analysis_col,first)
    return getattr(temp, second)
