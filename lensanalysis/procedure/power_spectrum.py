import numpy as np
import astropy.units as u

from lenstools.extern import _topology
from lenstools.utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

from .procedure import ConversionProcedureStep
from ..misc.log import logprocedure
from ..misc.feature_object import TomoPowerSpectra


class CalcTomoPowerSpectra(ConversionProcedureStep):
    """
    Computes the tomographic power spectra.

    Uses the same multipole bands for every bin combinations.

    Parameters
    ----------
    l_edges: array_like
        The edges of the multipole bands for every power spectrum. It must have
        shape (N,)
    scale: callable, optional
        Scaling to apply to the square of the Fourier pixels before harmonic 
        azimuthal averaging. Must be a function that takes the array of 
        multipole magnitudes as an input and returns an array of real numbers.
        Default is None. 

    Notes
    -----
    The implementation of conversion_operation is basically copied directly 
    from the cross bound method of the Spin0 class defined in 
    lenstools.image.convergence - this was done instead of simply calling the 
    method to reduce the number of fourier transforms by a factor of 
    (len(data_object)+1)).
    """

    def __init__(self,l_edges,scale=None):
        l_edges = np.asanyarray(l_edges)

        if len(l_edges.shape) != 1:
            raise ValueError("l_edges must have shape (N,)")
        elif l_edges.shape[0]<2:
            raise ValueError("l_edges must have at least 2 entries")
        elif (np.diff(l_edges)<=0):
            raise ValueError("l_edges must monotonically increase.")

        self.l_edges = l_edges

        if scale is None or callable(scale):
            self.scale = scale
        else:
            raise ValueError("scale must either be None or callable.")

    @classmethod
    def from_config(cls,ps_config):
        return cls(ps_config.get_tomo_multipole_bands())

    def conversion_operation(self,data_object,packet):

        scale = self.scale
        
        ft_maps = []
        sc_l = []
        side_angle = None
        shape = None
        for conv in data_object:

            # we probably do not need to perform the following error checking
            # if we are just piping the results of one procedure step to the
            # next
            if conv._masked:
                raise ValueError("Power spectrum calculation unable to handle "
                                 "masked maps.")
            elif conv.side_angle.unit.physical_type == "length":
                raise ValueError("Power spectrum calculation unable to handle "
                                 "when the side physical unit is 'length'")

            if side_angle is None:
                side_angle = conv.side_angle.to(u.deg).value
                shape = conv.data.shape
            else:
                if side_angle != conv.side_angle.to(u.deg).value:
                    raise ValueError("side_angles of convergence maps of "
                                     "different tomographic bins are not "
                                     "constant.")
                elif shape != conv.shape:
                    raise ValueError("the shapes of convergence maps of "
                                     "different tomographic bins are not "
                                     "constant.")

            ft_maps.append(fftengine.rfft2(conv.data))
            if scale is None:
                sc_l.append(None)
            else:
                sc_l.append(scale(conv.getEll()))

        l_edges = self.l_edges
        spectra = []

        for i,(ft_map1,sc) in enumerate(zip(ft_maps,sc_l)):
            for ft_map2 in ft_maps[i:]:
                spectra.appent(_topology.rfft2_azimuthal(ft_map1, ft_map2,
                                                         side_angle, l_edges,
                                                         sc))
        return [TomoPowerSpectra(spectra,len(data_object))]
