import warnings

import numpy as np
import astropy.table as tbl
import astropy.units as u
from lenstools import ShearMap

from .procedure import ConversionProcedureStep
from ..misc.log import logprocedure
from ..masked_operations import convert_shear_to_convergence, \
    determine_kernel_fft, convert_shear_to_smoothed_convergence_main

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
        Whether or not the resulting convergence map should be masked. Default 
        is False.
    pre_KS_smoothing: bool, optional
        Whether or not the smoothing should be performed in real space on the 
        Shear Map before performing the Kaiser-Squires Transorm. Default is 
        False.
    """

    def __init__(self,npixel,edge_angle,scale_angle,mask_result = False,
                 pre_KS_smoothing = False):
        assert npixel>0
        assert (edge_angle.unit.physical_type == "angle" and
                edge_angle.value > 0)
        assert (scale_angle.unit.physical_type == "angle" and
                scale_angle.value > 0)

        sigma_pix = (scale_angle * float(npixel)
                     / (edge_angle)).decompose().value
        fft_kernel, pad_axis = determine_kernel_fft(sigma_pix, npixel)
        self.fft_kernel = fft_kernel
        self.pad_axis = pad_axis
        self.npixel = npixel
        self.side_angle = edge_angle

        self.mask_result = mask_result
        self.pre_KS_smoothing = pre_KS_smoothing

        if self.mask_result and self.pre_KS_smoothing:
            raise NotImplementedError("Not currently equipped to mask result "
                                      "if we smoothing ahead of time.") 

    def conversion_operation(self,data_object,packet):
        logprocedure.debug(("Converting shear map(s) into smoothed convergence "
                            "map(s) for realization "
                            "{:d}").format(packet.data_id))

        out = []
        for shear_map in data_object:
            assert shear_map.side_angle == self.side_angle
            assert shear_map.data.shape[-1] == self.npixel
            conv = convert_shear_to_smoothed_convergence_main(shear_map,
                                                              self.pad_axis,
                                                              self.fft_kernel,
                                                              None, 0,
                                                              self.mask_result)
            out.append(conv)
        return out

class ConvMapToShearMap(ConversionProcedureStep):
    """
    Converts collections of convergence maps into collections of shear maps.
    """
    def conversion_operation(self,data_object,packet):
        out = [ShearMap.fromConvergence(conv_map) for conv_map in data_object]
        return out
