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

def _generate_seed(start_seed,map_id):
    if start_seed is None:
        return start_seed
    seed = start_seed + map_id
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

    May want to revist

    Parameters
    ----------
    start_seed : int
        Starting seed to which map_id is added. If this is None, then the seed 
        is not set. It needs to lie within the region of allowed seeds (the 
        range of a 32bit unsigned integer). Default is None.
    """

    def __init__(self,start_seed=None, inplace = False):
        _check_start_seed(start_seed)
        self.start_seed = start_seed
        self.inplace = False

    def add_noise(self,data_object,map_id):
        logprocedure.debug(("Adding noise to Shear Catalog(s) of realization "
                            "{:d}").format(map_id))
        seed = _generate_seed(self.start_seed,map_id)
        #return [catalog.shapeNoise(seed = seed) for catalog in data_object]
        warnings.warn("NOT ADDING NOISE PROPERLY - CURRENTLY PASSING THROUGH "
                      "THE NOISELESS OBJECT", RuntimeWarning)
        return data_object

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

class CatalogSourceEllipticityAdder(NoiseAdder):
    """
    Generate a collection of shear catalogs with intrinsic source ellipticities.
    """


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
