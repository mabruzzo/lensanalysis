class CosmologyCollection(object):
    """
    Contains collection of cosmology data.

    It's probably smart to lazily initialize individual cosmology collections.
    We may also want to cache cosmology collections.
    """
    def get_cosmology_data(self,name):
        pass

class SamplingCosmologyCollection(CosmologyCollection):
    """
    Collections of all of the sampled cosmologies.

    Much more standardized than the fiducial.
    """
    
    def streamlined_get_cosmology_data(self,name):
        """
        The function returns the cosmology data that uses the default source 
        configuration.
        """
        pass

class FiducialCosmologyCollection(CosmologyCollection):
    """
    May not need to worry about lazy initialization.
    """

    def get_cosmology_data(self,name):
        pass

    def get_fiducial_name(self):
        pass

    def get_fiducial_data(self):
        """
        Returns the data of the fiducial cosmology.

        This is just a convenience method that calls get_cosmology_data.
        """
        name = self.get_fiducial_name()
        return self.get_cosmology_data(name)




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


class MapProductCollection(object):
    def __init__(self):
        self.noiseless_map = None
        self.noisy_map = None
        self.smoothed_map = None
        self.smoothed_noisy_map = None


class FeatureProductCollection(object):
    """
    If we add more features, this will be extended.
    """
    def __init__(self):
        self.peak_locations = None
        self.peak_counts = None
        self.tomo_peak_locations = None
        self.tomo_peak_counts = None
