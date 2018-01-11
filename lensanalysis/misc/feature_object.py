class FeatureObject(object):
    """
    Represents features extracted from (mock) lensing surveys.
    """

class PeakLocations(FeatureObject):
    """
    Represents the convergence map peaks. Specifically includes locations, and 
    peak heights.

    Parameters
    ----------
    heights : array_like
        The heights of the peaks. Has shape (N,)
    locations : Quantity
        The locations of the peaks. Has shape (N,2) or (2,N)
    
    """
    def __init__(self, heights, locations):
        heights = np.array(heights)
        locations = np.array(locations)
        if len(heights.shape)!=1:
            raise ValueError("heights must be 1D")
        if len(locations.shape)!=2:
            raise ValueError("locations must be 2D")

        if locations.shape[0] == 2:
            locations = locations.T

        if locations.shape[0] != heights.shape[0]:
            raise ValueError("Locations must have a shape of (%d,2) or (2,%d)"
                             % (heights.shape[0], heights.shape[0]))

        self._heights = heights
        self._locations = locations

    @property
    def heights(self):
        return np.copy(self._heights)

    @property
    def locations(self):
        return np.copy(self._locations)

    def histogram(self, bins = 10, range = None, normed = False,
                  weigths = None, density = None):
        """
        Comput the histogram of the peak heights.

        Returns
        -------
        hist : PeakCounts
             An instance of PeakCounts that represents the histogram.
        bin_edges : array of dtype float
             Return the bin edges ``(length(hist)+1)``.

        See numpy.histogram for documentation on the Parameters and Notes.
        """
        hist, bin_edges = np.hist(self.heights, bins = bins, range = range,
                                  normed = normed, weights = weights,
                                  density = density)
        return PeakCounts(hist),bin_edges

    def __repr__(self):
        return 'PeakLocations(%s,\n              %s)' % (self._heights,
                                                         self._locations)

class PeakCounts(FeatureObject):
    """
    Represents the convergence map peak counts. 

    This class just tracks the number of counts per bin. It doesn't track 
    details about the bins.

    Parameters
    ----------
    counts : array_like
        The number of peaks per bin
    """

    def __init__(self,counts):
        self._counts = np.array(counts)

    @property
    def counts(self):
        return np.copy(self._counts)

    def __repr__(self):
        return 'PeakCounts(%s)' % self._counts.__repr__()
