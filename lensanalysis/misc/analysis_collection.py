
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
        self.tomo_feature_products = None


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
