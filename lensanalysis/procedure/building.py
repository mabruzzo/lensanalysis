import copy
from itertools import chain
import warnings

from .conversions import ShearCatalogToShearMap, ShearMapToSmoothedConvMap
from .noise_addition import NoiseAdditionStep, CatalogShapeNoiseAdder
from .peak_counting import LocatePeaks, BinPeaks
from .procedure import CompositeProcedureStep
from .io_step import SaveCollectionEntries
from .smoothing import ConvMapSmoothing
from ..misc.analysis_collection import get_analysis_col_value
from ..misc.enum_definitions import DescriptorEnum, all_combinations

def _equal_analysis_object(object_tuple1, object_tuple2):
    assert len(object_tuple1) == 2
    assert len(object_tuple2) == 2

    if object_tuple1[1] == object_tuple2[1]:
        if len(set(object_tuple1[0]) ^ set(object_tuple2[0])) == 0:
            return True
    return False


shear_cat = {}

_noisy = (DescriptorEnum.noisy,)
_smooth = (DescriptorEnum.smoothed,)
_noisy_smooth = (DescriptorEnum.noisy, DescriptorEnum.smoothed)

class AnalysisObjectList(object):
    def __init__(self,iterable = []):
        self._l = []
        for i in iterable:
            self.append(i)

    def __getitem__(self,index):
        return self._l[index]

    def __setitem__(self,index, value):
        self._l[index] = value

    def __delitem__(self,index):
        del self._l[index]

    def append(self, value):
        self._l.append(value)

    def remove(self,value):
        for elem in self._l:
            if _equal_analysis_object(elem, value):
                self._l.remove(elem)
                break
        else:
            raise ValueError(("{:s} is not contained within "
                              "AnalysisObjectList").format(str(value)))

    def __contains__(self,value):
        for elem in self._l:
            if _equal_analysis_object(elem, value):
                return True
        return False

    def __str__(self):
        return self._l.__str__()

    def __repr__(self):
        return "AnalysisObjectList({:s})".format(repr(self._l))

    def tolist(self):
        return copy.copy(self._l)

    def __len__(self):
        return len(self._l)
        
def _determine_save_analysis_objects(save_config):
    iterator = all_combinations(omit_analysis_objects = ["shear_cat",
                                                         "peak_counts", 
                                                         "peak_loc"],
                                omit_descriptors = [DescriptorEnum.tomo])
    iterator = chain(iterator,(((),"peak_counts"), ((),"peak_loc")))
    func = lambda x: get_analysis_col_value(save_config,*x)
    return AnalysisObjectList(filter(func,iterator))
    
def _check_build_peak_counting(save_config):
    """
    Check whether or not we should build anything related to peak counting.
    """
    if save_config.feature_products.peak_locations:
        return True
    elif save_config.feature_products.peak_counts:
        return True
    return False

def _wrap_save_step(decorator_step, save_step, following_sequential_step):
    """
    Helper function that wraps the save step and possibly the following 
    sequential step.
    """
    if following_sequential_step is None:
        decorator_step.wrapped_step = save_step
    else:
        temp = CompositeProcedureStep()
        temp.add(save_step)
        temp.add(following_sequential_step)
        decorator_step.wrapped_step = temp
    return decorator_step

def build_peak_counting(begin, procedure_config, save_config,
                        storage_collection, objects_to_save):
    if _equal_analysis_object(begin,((),"peak_counts")):
        raise ValueError("There is nothing to be done!")
    elif save_config.feature_products.peak_counts:
        # get the peak count bins
        peak_count_bins = procedure_config.peak_count_bins()
        # build the binning step
        second_step = BinPeaks(peak_count_bins)

        # check to see if we want to save
        if save_config.feature_products.peak_counts:
            storage = storage_collection.feature_products.peak_counts
            temp = SaveCollectionEntries(storage)
            second_step.wrapped_step = temp

            objects_to_save.remove(((),"peak_counts"))

        if _equal_analysis_object(begin,((),"peak_loc")):
            return second_step
    else:
        second_step = None    

    step = LocatePeaks()

    if save_config.feature_products.peak_locations:
        storage = storage_collection.feature_products.peak_locations
        save_step = SaveCollectionEntries(storage)
        objects_to_save.remove(((),"peak_loc"))
        _wrap_save_step(step, save_step, second_step)
    else:
        step.wrapped_step = second_step
    return step


def build_smooth_noisy_convergence(begin, procedure_config, save_config,
                                   storage_collection, objects_to_save,
                                   feature_step = None):
    """
    It will probably be better (faster) to do the convolution in Fouier Space.

    The only issue is that convolution theorem for discrete Fourier Transform 
    appears to be for circular convolution. I need to chat with Jose about it.
    """

    edge_angle = procedure_config.get_conv_map_edge_angle()
    npixel = procedure_config.get_conv_map_resolution()

    out = ShearMapToSmoothedConvMap(npixel,edge_angle,
                                    procedure_config.get_smoothing_scale())

    if save_config.conv_map.smoothed_noisy_map:
        storage = storage_collection.conv_map.smoothed_noisy_map
        objects_to_save.remove(((_noisy_smooth),"conv_map"))
        save_step = SaveCollectionEntries(storage)
        _wrap_save_step(out, save_step, feature_step)
    else:
        out.wrapped_step = feature_step
    return out

def build_shear_conversion(begin, procedure_config, save_config,
                           storage_collection, objects_to_save,
                           noisy = True, produce_single = True,
                           following_sequential_step = None):
    """
    Build the steps involving the conversion of the noiseless shear catalog to 
    a noisy (noiseless) shear map.
    """

    map_size = procedure_config.get_conv_map_edge_angle()
    npixel = procedure_config.get_conv_map_resolution()
    
    second_step = ShearCatalogToShearMap(map_size = map_size, npixel = npixel,
                                         produce_single = produce_single)

    if save_config.shear_map.noisy_map:
        # should probably factor this out as a function
        storage = storage_collection.shear_map.noisy_map
        objects_to_save.remove(((_noisy),"shear_map"))
        save_step = SaveCollectionEntries(storage)
        _wrap_save_step(second_step, save_step, following_sequential_step)
    else:
        second_step.wrapped_step = following_sequential_step

    # now we focus on whether or not noise was added
    if not noisy:
        return second_step

    warnings.warn("Currently the noise additon occurs forcing different seeds"
                  "for shear catalogs in the same realization, does not add "
                  "noise in place and does not include the rs_correction.",
                  RuntimeWarning)
    start_seed = procedure_config.get_noise_seed()
    out = NoiseAdditionStep(CatalogShapeNoiseAdder(start_seed,
                                                   diff_seed_elem = True,
                                                   rs_correction = False,
                                                   inplace = False))
    out.wrapped_step = second_step
    return out

        

def _build_procedure_helper(begin_object, procedure_config, save_config,
                            storage_collection, objects_to_save):

    # first we see if we need to build peak counting
    if _check_build_peak_counting(save_config):
        feature_step = build_peak_counting(begin_object, 
                                           procedure_config, save_config,
                                           storage_collection, 
                                           objects_to_save)
    else:
        feature_step = None

    if _equal_analysis_object(begin_object, (_noisy_smooth,"conv_map")):
        return feature_step

    # create the step to generate the smoothed noisy convergence map

    recent_step = build_smooth_noisy_convergence(begin_object, 
                                                 procedure_config,
                                                 save_config,
                                                 storage_collection,
                                                 objects_to_save,
                                                 feature_step)
    
    if _equal_analysis_object(begin_object, (_noisy, "shear_map")):
        return recent_step

    # create the step to generate the noisy shear map
    recent_step = build_shear_conversion(begin_object, procedure_config, 
                                         save_config,
                                         storage_collection, 
                                         objects_to_save,
                                         True, True, recent_step)

    if _equal_analysis_object(begin_object, ((),"shear_cat")):
        return recent_step
    else:
        raise ValueError("Not applicable beginning step")
    
def build_procedure(begin_object, procedure_config, save_config,
                    storage_collection):
    objects_to_save = _determine_save_analysis_objects(save_config)
    step = _build_procedure_helper(begin_object, procedure_config, save_config,
                                   storage_collection, objects_to_save)

    if len(objects_to_save)>0:
        raise ValueError(("Not all specified analysis products are scheduled "
                          "for saving. The remaining products to be saved are: "
                          "\n{:s}").format(str(objects_to_save)))
    return step
