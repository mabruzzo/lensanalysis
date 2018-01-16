from enum_definitions import DescriptorEnum

class AnalysisProductCollection(object):
    """
    Keeps track of the analysis products.
    """

    def __init__(self):
        self.conv_map = None
        self.tomo_conv_map = None
        self.shear_map = None
        self.tomo_shear_map = None
        self.feature_products = None


class ConvergenceMapProductCollection(object):
    def __init__(self):
        self.noiseless_map = None
        self.noisy_map = None
        self.smoothed_map = None
        self.smoothed_noisy_map = None

class ShearMapProductCollection(object):
    def __init__(self):
        self.noiseless_map = None
        self.noisy_map = None


class FeatureProductCollection(object):
    """
    If we add more features, this will be extended.
    """
    def __init__(self):
        self.peak_locations = None
        self.peak_counts = None
        self.tomo_peak_locations = None
        self.tomo_peak_counts = None



def _get_attribute_tuple(descriptors, object_name):
    if object_name == "conv_map":
        if DescriptorEnum.tomo in descriptors:
            first = 'tomo_conv_map'
        else:
            first = 'conv_map'

        if DescriptorEnum.noisy in descriptors:
            if DescriptorEnum.smooth in descriptors:
                second = 'smoothed_noisy_map'
            else:
                second = 'noisy_map'
        elif DescriptorEnum.smooth in descriptors:
            second = 'smoothed_map'
        else:
            second = 'noiseless_map'

    elif object_name == "shear_map":
        if DescriptorEnum.tomo in descriptors:
            first = 'tomo_shear_map'
        else:
            first = 'shear_map'

        if DescriptorEnum.noisy in descriptors:
            second = 'noisy_map'
        else:
            second = 'noiseless_map'

    elif object_name in ["peak_loc","peak_counts"]:
        first = 'feature_products'

        if object_name == "peak_loc":
            if DescriptorEnum.tomo in descriptors:
                second = "tomo_peak_locations"
            else:
                second = "peak_locations"
        if DescriptorEnum.tomo in descriptors:
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
