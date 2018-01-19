from procedure import ConversionProcedureStep
from ..misc.log import logprocedure

class LocatePeaks(ConversionProcedureStep):
    """
    Locates the peaks of the convergence map.
    """

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Locating peaks for realization "
                            "{:d}").format(packet.data_id))
        out = []
        for elem in data_object:
            heights,positions = elem.locatePeaks(extent)
            out.append(PeakLocations(heights = heights, locations = positions))
        return out

class BinPeaks(ConversionProcedureStep):
    """
    Bins the peaks of the convergence map.
    
    Configureable to work with tomography or a single grouping. (will need to 
    revisit)

    Parameters
    ----------
    bin_edges : array_like
        The edges of the bins for the peak count histogram.
    """
    def __init__(self,bin_edges):
        bin_edges = np.array(bin_edges)
        if len(bin_edges.shape)!=1:
            raise ValueError("bin_edges must be 1 dimensional")
        self.bin_edges = bin_edges

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Binning peaks for realization "
                            "{:d}").format(packet.data_id))
        if len(data_object)>1:
            raise RuntimeError("Need to update to work with different bins")
        out = []
        for peak_loc in data_object:
            hist, bin_edges = peak_loc.histogram(self.bin_edges)
            out.append(hist)
        return out
