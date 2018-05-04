"""
Possibly Useful Notes:
    There appear to be 2 main types of features with respect to tomography:
    1. Generalizable - The tomographic version is simply a list of the features
       for individual tomographic bins. Examples include the Peaks and moments
    2. CrossStatistic - The tomographic version of a feature includes 
       quantification of characteritics relative to the characteristics in 
       other bins in addition to the quanitification of the characteristic in 
       its own bin. Therefore the representation of the Tomographic and 
       non-Tomographic versions are different. An example of this is the Power 
       Spectrum.

I think it still makes sense to store feature objects in iterable containers as 
we have previously done. For now, we will store CrossStatistics (whether they 
are for tomography or not) in lists while passing them around through the 
procedure object, but it may be worthwhile to juse pass around these 
CrossStatistic objects themselves and give them the interface of the containers.
"""

import numpy as np
import astropy.units as u

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
        assert isinstance(locations,u.Quantity)
        locations = u.Quantity(locations)
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
        return self._heights.copy()

    @property
    def locations(self):
        return self._locations.copy()

    def histogram(self, bins = 10, range = None, normed = False,
                  weights = None, density = None):
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
        hist, bin_edges = np.histogram(self.heights, bins = bins, 
                                       range = range, normed = normed, 
                                       weights = weights, density = density)
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


class TomoPowerSpectra(FeatureObject):
    """
    Represents the Tomographic Convergence Power Spectra

    Parameters
    ----------
    power : array_like
        The power for each multipole band in every power spectra. Has shape 
        (N,M) where N is the number of spectra, given by 
        (num_tomo_bins+1)*num_tomo_bins//2, and M is the number of 
        multipole bands per spectrum.
    num_tomo_bins : int
        The number of tomographic bins.

    Notes
    -----
    Because the Tomographic Power Spectra is a cross statistic - the 
    multipole moments of the fourier transformed Convergence Maps are 
    correlated between all unique combinations of the tomographic bins - this 
    can not be simply viewed as a list of features from each individual bin.
 
    To further elaborate, in the case of a single bin, only one power 
    spectrum is computed per each bin. In the tomographic case, if there are 
    are n_tomo tomographic bins, then num_tomo_bins*(num_tomo_bins+1)/2 power 
    spectrum are computed (including the autocorrelated power spectra for each 
    bin).

    Internally, we keep track of the various Spectra in a scheme resembling 
    an upper triangular matrix with locations corresponding to the unique 
    combinations of various tomographic bins (using zero-indexing).

    let n = (num_tomo_bins-1)

    P(l)_ij = [[P(l)_0,0    P(l)_0,1    P(l)_0,2       ...      P(l)_0,n  ]
               [            P(l)_1,1    P(l)_1,2       ...      P(l)_1,n  ]
               [                           ...         ...         ...    ]
               [                                       ...      P(l)_n-1,n]
               [                                                P(l)_n,n  ]]

    """

    def __init__(self,power, num_tomo_bins):
        if not isinstance(num_tomo_bins,int):
            raise ValueError("num_tomo_bins must be an int.")
        elif num_tomo_bins <1:
            raise ValueError("num_tomo_bins must be a positive integer")

        self.num_tomo_bins = num_tomo_bins
        power = np.asanyarray(power)
        if len(power.shape)!=2:
            raise ValueError("power must be a 2D array")
        n,m = power.shape

        expected = self.num_tomo_bins*(self.num_tomo_bins+1)//2
        if n!= expected:
            raise ValueError(("Axis 0 of power must have a length of {:d} "
                              "if num_tomo_bins = {:d}").format(expected,
                                                                num_tomo_bins))
        if m == 0:
            raise ValueError("Axis 1 of power must have non-zero length")
        self._power = power

    def get_spectrum(self,i,j):
        """
        Returns the spectrum from crossing bin i with bin j (zero-indexed).

        Notes
        -----
        We take advantage of the fact that the number of spectra is a triangle 
        number: index = (self.num_tomo_bins *(self.num_tomo_bins+1)//2 
                         - (self.num_tomo_bins-i) *(self.num_tomo_bins-i+1) //2 
                         + j - i)
        """
        if i < 0 or i >= self.num_tomo_bins:
            raise ValueError("i must be at least 0 but less than "
                             "{:d}".format(self.num_tomo_bins))
        if j < i or j >= self.num_tomo_bins:
            raise ValueError("j must be at least {:d} but less than ".format(i)
                             "{:d}".format(self.num_tomo_bins))
        return np.copy(self._power[i*(2*self.num_tomo_bins + 1 - i)//2 + j - i,
                                   :])

    @property
    def power(self):
        return np.copy(self._power)
    
    def get_all_features(self):
        return self.power.flatten(order='C')

    def __repr__(self):
        return ('TomoPowerSpectrum({:s},\n'
                '                  {:d})').format(self._power.__repr__(),
                                                  self.num_tomo_bins)
