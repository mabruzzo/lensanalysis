"""

Defining classes to add noise.

 - Might want to create look up tables in which we can register additional 
   noise adders.
 - May also want to give additional thought to how we generate seeds 
 - May want to add a NoiseAdder that loads previously computed noise maps

"""

from abc import ABCMeta, abstractmethod
import warnings

import numpy as np
import astropy.table as tbl

from procedure import IntermediateProcedureStep
from ..misc.log import logprocedure

class NoiseAdder(object):
    __metaclass__ = ABCMeta

    """
    Abstract Base Class for objects that add noise to ShearCatalogs
    """

    @abstractmethod
    def add_noise(self,data_object,map_id):
        pass

_uii32 = np.iinfo(np.uint32)

def _check_start_seed(start_seed):
    if start_seed is not None:
        if not isinstance(start_seed,int):
            raise TypeError("start_seed must be an int")
        if start_seed<0:
            raise ValueError("start_seed must be greater than or equal to "
                             "0.")
        elif start_seed>_uii32.max:
            raise ValueError(("start_seed can be no larger that %d,\nthe "
                              "maximum value of a 32 bit unsigned integer.")
                             % int(_uii32.max))

def _generate_seed(start_seed,map_id,elem_id,num_elem):
    if start_seed is None:
        return start_seed
    if num_elem == 0:
        delta = map_id
    else:
        delta = map_id*num_elem + elem_id

    seed = start_seed + delta
    if seed > _uii32.max:
        seed = (seed%_uii32.max)
    return seed

class CatalogShapeNoiseAdder(NoiseAdder):
    """
    Generates a catalog collection with randomly drawn shape noise.

    Presently we start with a starting seed and just add the number of the 
    realization to it. If the sum of the numbers would be greater than the 
    maximum value of an 32bit unsigned integer (the maximum value in the range 
    of seeds, we have it wrap around)

    A single instance of this object is only meant to be used to add noise to 
    collections of shear maps of fixed size. If we choose to have unique 
    seeds for each shear map, then we internally track the size of the 
    collections the first time we actually add noise. An error will be raised 
    if any further collection has a different number of items.

    Parameters
    ----------
    start_seed : int
        Starting seed to which map_id is added. If this is None, then the seed 
        is not set. It needs to lie within the region of allowed seeds (the 
        range of a 32bit unsigned integer). Default is None.
    diff_seed_elem : bool
        Whether or not to use a different seed for different shear maps in a 
        data_collection. If this is true, then we strictly require that all 
        catalog collections we add noise to are the same size. Default is 
        False.
    rs_correction : bool
        Include denominator (1+g*es) in the shear correction. Default is True.
    inplace : bool
        If the data objects should be modified in place. Default is False.
    """

    def __init__(self,start_seed=None, diff_seed_elem = False, 
                 rs_correction = True, inplace = False):
        _check_start_seed(start_seed)
        self.start_seed = start_seed
        self.inplace = inplace
        
        self.diff_seed_elem = diff_seed_elem
        if diff_seed_elem:
            self._num_elem = None
        else:
            self._num_elem = 0

        self.rs_correction = rs_correction

    def add_noise(self,data_object,map_id):
        logprocedure.debug(("Adding noise to Shear Catalog(s) of realization "
                            "{:d}").format(map_id))
        start_seed = self.start_seed
        inplace = self.inplace
        if self._num_elem is None:
            self._num_elem = len(data_object)
        elif self._num_elem>0:
            assert number_elem == len(data_object)
        number_elem = self._num_elem
        rs_correction = self.rs_correction


        out = []
        for i,catalog in enumerate(data_object):
            seed = _generate_seed(start_seed,map_id,i,number_elem)
            noise = catalog.shapeNoise(seed)
            if not inplace:
                # we are trying to trying to maintain the original data 
                # without making copies
                result = catalog.copy(copy_data = False)
                result.remove_columns(["shear1","shear2"])
            else:
                result = catalog

            temp =catalog.addSourceEllipticity(noise, 
                                                 es_colnames = ("shear1",
                                                                "shear2"),
                                                 rs_correction = rs_correction,
                                                 inplace = False)
            result["shear1"] = temp["shear1"]
            result["shear2"] = temp["shear2"]
            out.append(result)
        return out

class ConvMapNormalShapeNoiseAdder(NoiseAdder):
    """
    Generates a catalog collection with randomly drawn shape noise from a 
    normal distribution.

    Presently we start with a starting seed and just add the number of the 
    realization to it. If the sum of the numbers would be greater than the 
    maximum value of an 32bit unsigned integer (the maximum value in the range 
    of seeds, we have it wrap around)

    May want to revist - (For example if performing tomography with 
    convergence maps at fixed locations).

    Parameters
    ----------
    mean : float
        The mean of the normal distribution. Default is 0.0.
    scale : float
        The standard deviation of the normal distribution.
    start_seed : int
        Starting seed to which map_id is added. If this is None, then the seed 
        is not set. It needs to lie within the region of allowed seeds (the 
        range of a 32bit unsigned integer). Default is None.
    """
    def __init__(self,mean=0.0,scale=1.0,start_seed=None):
        _check_start_seed(start_seed)
        self.mean = mean
        self.scale = scale
        self.start_seed = start_seed

    def add_noise(self,data_object,map_id):
        seed = _generate_seed(start_seed,map_id)
        if seed is not None:
            np.random.seed(seed)
        out = []
        for conv_map in data_object:
            shape = conv_map.data.shape
            noise = np.random.normal(loc = self.mean, scale = self.scale,
                                     shape = shape)
            out.append(conv_map + noise)
        return out


class NoiseAdditionStep(IntermediateProcedureStep):
    """
    Adds shape noise to the entries of the data object collection.
    """

    def __init__(self,noise_adder=None):
        self.noise_adder = noise_adder

    def intermediate_operation(self,data_object,packet):
        
        if self.noise_adder is None:
            return data_object
        else:
            return self.noise_adder.add_noise(data_object,packet.data_id)
