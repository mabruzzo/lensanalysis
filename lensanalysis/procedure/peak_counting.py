from numbers import Number

import numpy as np

from .procedure import ConversionProcedureStep
from ..misc.log import logprocedure
from ..misc.feature_object import PeakLocations

class LocatePeaks(ConversionProcedureStep):
    """
    Locates the peaks of the convergence map.

    Parameters
    ----------
    norm : boolean
        Indicates if the peak counts should be normalized by the standard 
        deviation of the map (effectively causes the peak heights to be SNR 
        where sigma is unique to each map).
    """

    def __init__(self,norm = False):
        self.norm = norm

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Locating peaks for realization "
                            "{:d}").format(packet.data_id))
        out = []
        for elem in data_object:
            extent = [np.nanmin(elem.data), np.nanmax(elem.data)]
            heights,positions = elem.locatePeaks(extent,norm=self.norm)
            out.append(PeakLocations(heights = heights, locations = positions))
        return out

def _contains_1D_arrays(sequence):
    out = []
    for elem in sequence:
        array = np.array(elem)
        if len(array.shape)!=1:
            raise ValueError("each array_like object must be 1 dimensional")
        out.append(array)
    return out

class BinPeaks(ConversionProcedureStep):
    """
    Bins the peaks of the convergence map.
    
    Configureable to work with tomography or a single grouping. (will need to 
    revisit)

    Parameters
    ----------
    bin_edges : array_like or Sequence of array_like objects
        The edges of the bins for the peak count histogram. If this is to work 
        with tomography, it must a sequences of the edges of the bins for each  
        peak histogram.
    """
    def __init__(self,bin_edges):
        if isinstance(bin_edges[0],Number):
            bin_edges = np.array(bin_edges)
            if len(bin_edges.shape)!=1:
                raise ValueError("bin_edges must be 1 dimensional")
            self.num_tomo_bins = 0
        else:
            bin_edges = _contains_1D_arrays(bin_edges)
            if len(bin_edges) == 0:
                raise ValueError("need to include at least 1 bin")
            elif len(bin_edges) == 1:
                bin_edges = bin_edges[0]
                self.num_tomo_bins = 0
            else:
                self.num_tomo_bins = len(bin_edges)
        self.bin_edges = bin_edges

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Binning peaks for realization "
                            "{:d}").format(packet.data_id))
        out = []
        if self.num_tomo_bins:
            if len(data_object) != self.num_tomo_bins:
                raise ValueError(("BinPeaks was configured to work with {:d}"
                                  "tomographic bins. Instead,\n{:d} bins are "
                                  "being processed").format(self.num_tomo_bins,
                                                            len(data_object)))
            for peak_loc, bin_edges in zip(data_object,self.bin_edges):
                hist, bin_edges = peak_loc.histogram(bin_edges)
                out.append(hist)
        else:
            if len(data_object)>1:
                raise ValueError("BinPeaks was not configured to work with "
                                 "tomographic data")
            for peak_loc in data_object:
                hist, bin_edges = peak_loc.histogram(self.bin_edges)
                out.append(hist)
        return out
