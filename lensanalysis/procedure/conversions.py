import astropy.table as tbl
from lenstools import ShearMap

from procedure import ConversionProcedureStep

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

        out = []
        if self.produce_single:
            catalogs = [tbl.vstack(data_object)]
        else:
            catalogs = data_object

        for catalog in catalogs:
            out.append(catalog.toMap(mapsize = self.map_size,
                                     npixel = self.npixel, smooth = self.smooth,
                                     **self.kwargs))
        return out

class ShearMapToConvMap(ConversionProcedureStep):
    """
    Converts collections of shear maps into collections of convergence maps.
    """
    def conversion_operation(self,data_object,packet):
        out = [shear_map.convergence() for shear_map in data_object]
        return out

class ConvMapToShearMap(ConversionProcedureStep):
    """
    Converts collections of convergence maps into collections of shear maps.
    """
    def conversion_operation(self,data_object,packet):
        out = [ShearMap.fromConvergence(conv_map) for conv_map in data_object]
        return out
