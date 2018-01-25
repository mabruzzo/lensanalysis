from procedure import IntermediateProcedureStep
from ..misc.log import logprocedure
from ..masked_operations import smooth_conv_map

class ConvMapSmoothing(IntermediateProcedureStep):
    """
    Produces a collection of smoothed convergence maps.

    Parameters
    ----------
    scale_angle : astropy.Quantity
        Size of smoothing to be performed. Must have units.
    kind : str
        Type of smoothing to be performed. Options are "gaussian" and 
        "gaussianFFT". Default is "gaussian."
    kwargs : dict
        Keyword arguments to be passed to the filter function.
    """
    def __init__(self,scale_angle, kind = "gaussian", **kwargs):
        self.scale_angle = scale_angle
        self.kind = "gaussianfft"
        self.kwargs = kwargs

    def intermediate_operation(self,data_object,packet):
        logprocedure.debug("Smoothing realization {:d}".format(packet.data_id))
        out = []
        for conv_map in data_object:
            out.append(smooth_conv_map(conv_map, self.scale_angle,
                                       kind = self.kind, **self.kwargs))
        return out
