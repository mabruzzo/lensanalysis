import warnings

import numpy as np
import astropy.table as tbl
import astropy.units as u
from scipy.ndimage import filters
from lenstools import ShearMap

from .procedure import ConversionProcedureStep
from ..misc.log import logprocedure
from ..masked_operations import convert_shear_to_convergence, \
    determine_kernel_fft, smooth_shear, clip_smoothed_conv_map_boundaries, \
    determine_kernel_length
from ..masked_operations import convert_shear_to_smoothed_convergence_main \
    as convert_sh_to_sm_conv_main

class ShearCatalogToShearMap(ConversionProcedureStep):
    """
    Converts collections of shear catalogs into collections of shear maps.

    Parameters
    ----------
    map_size : astropy.Quantity
        Spatial size of each map
    npixel : int
        Number of pixels on a side of each map
    smooth : astropy.Quantity
        If not None, each map is smoothed with a gaussian filter of scale 
        smooth. Default is None.
    produce_single : bool
        Whether one shear map (True) should be produced from all shear 
        catalogs in a shear catalog collection or one shear map should be 
        produced for each shear catalog in the shear catalog collection. 
        Default is True.
   kwargs : dict
        Additional arguments to be passed to the pixelize method of each shear 
        catalog being smoothed.
    
    """
    def __init__(self, map_size, npixel, smooth = None, produce_single = True,
                 **kwargs):
        self.map_size = map_size
        self.npixel = npixel
        self.smooth = smooth
        self.produce_single = produce_single
        self.kwargs = kwargs

    def conversion_operation(self,data_object,packet):
        if self.produce_single:
            message = ("Pixelizing Shear Catalog(s) of realization {:d} into "
                       "a single shear map.")
        else:
            message = ("Pixelizing Shear Catalog(s) of realization {:d} into "
                       "shear maps.")
        logprocedure.debug(message.format(packet.data_id))

        out = []
        if self.produce_single:
            catalogs = [tbl.vstack(data_object)]
        else:
            catalogs = data_object

        for catalog in catalogs:
            out.append(catalog.toMap(map_size = self.map_size,
                                     npixel = self.npixel,
                                     smooth = self.smooth,
                                     **self.kwargs))
        return out

class ShearMapToConvMap(ConversionProcedureStep):
    """
    Converts collections of shear maps into collections of convergence maps.
    """
    def conversion_operation(self,data_object,packet):
        raise RuntimeError("Need to implement option about whether or not the "
                           "resulting map should be remasked")
        logprocedure.debug(("Converting shear map(s) into convergence map(s) "
                            "for realiztion {:d}").format(packet.data_id))

        out = map(convert_shear_to_convergence, data_object)
        return out


def _mask_from_shear(shear_map,conv_map):
    """ 
    Helper Function to apply shear_map mask on smoothed convergence map inplace.
    """
    m = np.logical_or(np.isnan(shear_map.data[0],
                               shear_map.data[1]))
    conv_map.mask(m.astype(np.int8),inplace=True)
    
class ShearMapToSmoothedConvMap(ConversionProcedureStep):
    """
    Converts collections of shear maps into collections of smoothed convergence 
    maps.

    Parameters
    ----------
    npixel : int
        The number of pixels along the edge of the convergence map
    edge_angle : astropy.Quantity
        The angle along the edge of the convergence map.
    scale_angle : astropy.Quantity
        Size of smoothing to be performed. Must have units.
    mask_result : bool, optional
        Whether or not the resulting convergence map should be masked. Note 
        that the original shear map is always unmasked and has masked values 
        set to 0. Only at the very end, after smoothing and applying the 
        Kaiser-Squires Transform is the convergence map masked. Default is
        False.
    real_space_smoothing : bool, optional
        If the map should be smoothed in real space. Default is False.
    edge_mode: str,optional
        How to handle smoothing at the edge of the array. Valid modes are 
        {'constant', 'mirror'}. When mode is equal to 'constant' it assumes 
        everywhere beyond the edge is set to 0. Presently, fft smoothing only 
        works with 'constant'. If set to None, then edge_mode defaults to 
        'mirror' for real space smoothing and 'constant' for Fourier space 
        smoothing.
    clip_boundaries: bool,optional
        If the boundaries of the convergence map where the smoothing kernel 
        required information about pixels found outside of the map should be 
        discarded. Default is False.

    Notes
    -----
        If smoothing is performed in Fourier Space, then the Kaiser-Squires 
        Transform is implicitly performed using an array that has been 
        zero-padded enough to avoid the wrapping boundary condition while 
        smoothing

        If smoothing is performed in Real Space, then there is no zero-padding
        for the Kaiser-Squires Transform.

        Additionally, the Kaiser-Squires Transform commutes with smoothing. 
        The exact same is achieved regardless of the order of operations. 
        Obviously, we will always smooth the convergence map.
    """

    def __init__(self,npixel,edge_angle,scale_angle,mask_result = False,
                 real_space_smoothing = False,
                 edge_mode = None, clip_boundaries = False):
        assert npixel>0
        assert (edge_angle.unit.physical_type == "angle" and
                edge_angle.value > 0)
        assert (scale_angle.unit.physical_type == "angle" and
                scale_angle.value > 0)

        sigma_pix = (scale_angle * float(npixel)
                     / (edge_angle)).decompose().value
        self.sigma_pix = sigma_pix

        if real_space_smoothing:
            if edge_mode is None:
                edge_mode = 'mirror'
            elif edge_mode not in ['mirror', 'constant']:
                raise ValueError("For real space smoothing, edge_mode must be "
                                 "'mirror' or 'constant'")
            
        else:
            if edge_mode is None:
                edge_mode = 'constant'
            elif edge_mode != 'constant':
                raise ValueError("For Fourier Space Smoothing, edge_mode must "
                                 "be 'constant'")
            fft_kernel, pad_axis = determine_kernel_fft(sigma_pix, npixel)
            self.fft_kernel = fft_kernel
            self.pad_axis = pad_axis

        self.npixel = npixel
        self.side_angle = edge_angle

        self.mask_result = mask_result
        self.real_space_smoothing = real_space_smoothing
        self.edge_mode = edge_mode

        if clip_boundaries:
            self.kernel_width = determine_kernel_length(sigma_pix,
                                                         truncate = 4.0)
            assert self.kernel_width > 0
        else:
            self.kernel_width = 0

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Converting shear map(s) into smoothed convergence "
                            "map(s) for realization "
                            "{:d}").format(packet.data_id))

        out = []
        real_space_smoothing = self.real_space_smoothing
        mode = self.edge_mode
        sigma_pix = self.sigma_pix
        kernel_width = self.kernel_width

        for shear_map in data_object:
            assert shear_map.side_angle == self.side_angle
            assert shear_map.data.shape[-1] == self.npixel
            if real_space_smoothing:
                # unmask shear map and convert to convergence map 
                temp = convert_shear_to_convergence(shear_map,
                                                    map_mask = None,
                                                    fill = 0,
                                                    perserve_mask = False)
                # now smooth the Convergence Map
                kwargs = dict((k,getattr(self,k)) \
                              for k in temp._extra_attributes)
                conv = ConvergenceMap(filters.gauusian_filter(temp.data,
                                                              sigma_pix,
                                                              mode = mode),
                                      temp.side_angle,**kwargs)
                # Finally, optionally mask the Convergence Map
                if self.mask_result:
                    _mask_from_shear(shear_map,conv)
            else:
                # as of now, we only allow smoothing in Fourier space using
                # the 'constant' edge_mode
                conv = convert_sh_to_sm_conv_main(shear_map, self.pad_axis,
                                                  self.fft_kernel, None, 0,
                                                  self.mask_result)
            if kernel_width:
                out.append(clip_smoothed_conv_map_boundaries(conv,
                                                             kernel_width))
            else:
                out.append(conv)
        return out

class ConvMapToShearMap(ConversionProcedureStep):
    """
    Converts collections of convergence maps into collections of shear maps.
    """
    def conversion_operation(self,data_object,packet):
        out = [ShearMap.fromConvergence(conv_map) for conv_map in data_object]
        return out
