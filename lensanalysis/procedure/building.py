import copy

from .conversions import ShearCatalogToShearMap, ShearMapToSmoothedConvMap
from .noise_addition import NoiseAdditionStep, CatalogShapeNoiseAdder
from .peak_counting import LocatePeaks, BinPeaks
from .procedure import CompositeProcedureStep
from .io_step import SaveCollectionEntries
from .smoothing import ConvMapSmoothing
from .tomo_binning import PseudoPhotozRebinner, ShearCatRebinning, \
    LazyStaticShearCatRebinning
from .power_spectrum import CalcTomoPowerSpectra
from ..misc.analysis_collection import default_value_UAPC
from ..misc.enum_definitions import Descriptor

def _convert_descriptor(descr_input):
    if isinstance(descr_input,Descriptor):
        return descr_input
    elif not isinstance(descr_input,(list,tuple)):
        raise TypeError

    descr = Descriptor.none
    for elem in descr_input:
        descr = descr | elem
    return descr


def _check_build_peak_counting(save_config,tomo=False):
    """
    Check whether or not we should build anything related to peak counting.
    """
    if tomo:
        return (save_config.feature_products.tomo_peak_locations or
                save_config.feature_products.tomo_peak_counts)
    else:
        return (save_config.feature_products.peak_locations or
                save_config.feature_products.peak_counts)

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

def build_peak_counting(begin, procedure_config, storage_collection,
                        objects_to_save, tomo = False, num_tomo_bins = 0):
    fp_storage = storage_collection.feature_products

    if tomo:
        save_peak_counts = objects_to_save.feature_products.tomo_peak_counts
        peak_count_storage = fp_storage.tomo_peak_counts
        save_peak_loc = objects_to_save.feature_products.tomo_peak_locations
        peak_loc_storage = fp_storage.tomo_peak_locations
        descriptor = Descriptor.tomo
    else:
        save_peak_counts = objects_to_save.feature_products.peak_counts
        peak_count_storage = fp_storage.peak_counts
        save_peak_loc = objects_to_save.feature_products.peak_locations
        peak_loc_storage = fp_storage.peak_locations
        descriptor = Descriptor.none

    if begin == (descriptor,"peak_counts"):
        raise ValueError("There is nothing to be done!")
    elif save_peak_counts:
        # get the peak count bins
        if tomo:
            peak_count_bins = procedure_config.tomo_peak_count_bins(\
                                                                num_tomo_bins)
        else:
            peak_count_bins = procedure_config.peak_count_bins()

        # build the binning step
        second_step = BinPeaks(peak_count_bins)

        # we can only want to save this
        if save_peak_counts:
            temp = SaveCollectionEntries(peak_count_storage)
            second_step.wrapped_step = temp
            if tomo:
                objects_to_save.feature_products.tomo_peak_counts = False
            else:
                objects_to_save.feature_products.peak_counts = False

        if begin == (descriptor,"peak_loc"):
            return second_step
    else:
        second_step = None

    normalize_sigma = procedure_config.sigma_indiv_map(tomo=tomo)
    step = LocatePeaks(norm=normalize_sigma)

    if save_peak_loc:
        save_step = SaveCollectionEntries(peak_loc_storage)
        if tomo:
            objects_to_save.feature_products.tomo_peak_locations = False
        else:
            objects_to_save.feature_products.peak_locations = False
        _wrap_save_step(step, save_step, second_step)
    else:
        step.wrapped_step = second_step
    return step

def build_power_spectrum(begin, procedure_config, storage_collection,
                         objects_to_save, tomo = False):
    if not tomo:
        raise NotImplementedError("We have not implemented non-tomographic "
                                  "power spectra.")
    try:
        config = procedure_config.feature_config['power_spectrum']
    except KeyError:
        raise ValueError("Power spectra section of procedure configuration is "
                         "not specified.")

    out = CalcTomoPowerSpectra.from_config(ps_config=config)

    if objects_to_save.feature_products.tomo_power_spectrum:
        storage = storage_collection.feature_products.tomo_power_spectrum
        save_step = SaveCollectionEntries(storage)
        out.wrapped_step = save_step
        objects_to_save.feature_products.tomo_power_spectrum = False
    else:
        # sanity check
        raise RuntimeError("There is no point in constructing the power "
                           "spectrum if we don't save it")
    return out

def build_smooth_noisy_convergence(begin, procedure_config, storage_collection,
                                   objects_to_save, feature_step = None,
                                   tomo = False):
    """
    Set's up conversion from shear map to smoothed convergence map
    """

    if tomo:
        ots_conv = objects_to_save.tomo_conv_map
        conv_map_storage = storage_collection.tomo_conv_map
    else:
        ots_conv = objects_to_save.conv_map
        conv_map_storage = storage_collection.conv_map

    if feature_step is None and not ots_conv.smoothed_noisy_map:
        return None

    edge_angle = procedure_config.get_conv_map_edge_angle()
    npixel = procedure_config.get_conv_map_resolution()
    mask_conv_map = procedure_config.mask_convergence_conversion()
    clip_boundaries = procedure_config.clip_boundaries()

    out = ShearMapToSmoothedConvMap(npixel,edge_angle,
                                    procedure_config.get_smoothing_scale(),
                                    mask_result = mask_conv_map,
                                    clip_boundaries = clip_boundaries)

    if ots_conv.smoothed_noisy_map:
        # get the storage object
        storage = conv_map_storage.smoothed_noisy_map

        # set the appropriate value of objects_to_save to False
        ots_conv.smoothed_noisy_map = False

        # build the actual procedural steps
        save_step = SaveCollectionEntries(storage)
        _wrap_save_step(out, save_step, feature_step)
    else:
        out.wrapped_step = feature_step
    return out

def build_shear_conversion(begin, procedure_config, storage_collection,
                           objects_to_save, following_sequential_step = None,
                           tomo =False):
    """
    Build the steps involving the conversion of the shear catalog to a shear 
    map.
    """
    produce_single = not tomo

    if tomo:
        ots_shear = objects_to_save.tomo_shear_map
        shear_map_storage = storage_collection.tomo_shear_map
    else:
        ots_shear = objects_to_save.shear_map
        shear_map_storage = storage_collection.shear_map

    if following_sequential_step is None and not ots_shear.noisy_map:
        return None

    map_size = procedure_config.get_conv_map_edge_angle()
    npixel = procedure_config.get_conv_map_resolution()

    real_space = procedure_config.real_space_smoothing()
    edge_mode = procedure_config.edge_mode()

    second_step = ShearCatalogToShearMap(map_size = map_size, npixel = npixel,
                                         produce_single = produce_single,
                                         real_space_smoothing = real_space,
                                         edge_mode = edge_mode)

    if ots_shear.noisy_map:
        # get the storage object
        storage = shear_map_storage.noisy_map

        # set the appropriate value of objects_to_save to False
        ots_shear.noisy_map = False

        # build the actual procedural steps
        save_step = SaveCollectionEntries(storage)
        _wrap_save_step(second_step, save_step, following_sequential_step)
    else:
        second_step.wrapped_step = following_sequential_step

    return second_step

def build_rebinning(begin, procedure_config, storage_collection,
                    objects_to_save, following_sequential_step = None,
                    num_tomo_bins=-1, ppz_noise = None):
    """
    Constructs the steps to rebin the Shear Catalog.
    """

    if procedure_config.has_non_ppz_rebinning():
        if ppz_noise is not None:
            raise RuntimeError("Not presently equipped to handle rebinning and "
                               "pseudo photoz errors simultaneously.")

        # first lets get the bin limits
        num_bins,rebinned_z_intervals = procedure_config.get_bin_limits()
        if rebinned_z_intervals is None:
            raise ValueError("Need the specific bin_intervals for rebinning.")

        # check for possible specification of a separate z catalog for
        # identifying mapping galaxies to the result tomography bin.
        temp = procedure_config.get_z_binning_cat()
        z_binning_cat_root_dir, z_binning_cat_fname_formatter = temp

        #for now we assume share_position_component is True
        share_position_component = True
        next_step = LazyStaticShearCatRebinning(rebinned_z_intervals,
                                                share_position_component,
                                                z_binning_cat_fname_formatter,
                                                z_binning_cat_root_dir)
        print next_step
    else:
        if ppz_noise is None:
            return following_sequential_step

        # first lets get values to build the rebinner - checks for the bin
        # limits
        temp = procedure_config.get_bin_limits()
        if temp is None or temp[0] is None:
            # in this case, the bin_intervals are provided
            num_bins,bin_intervals = procedure_config.get_bin_limits()
            if num_bins != num_tomo_bins:
                raise ValueError("For PseudoPhotoz using the original bins,\n"
                                 "the num_bins option of the Rebinning section "
                                 "of\nthe Procedure Configuration file must be "
                                 "-1\n(unspecified) or match the value of the\n"
                                 "num_tomo_bins option in the General section "
                                 "of\nthe storage configuration file.")
            kwargs = {}
        else:
            # in this case the bin_intervals are not provided. Therefore, we
            # need to check for guidelines to help determine the bin intervals
            bin_intervals = None
            kwargs = procedure_config.input_bin_lim_guidelines()
            # now let's construct the rebinner.
        rebinner = PseudoPhotozRebinner(ppz_noise,
                                        bin_intervals = bin_intervals,
                                        colname='z', **kwargs)
    
        next_step = ShearCatRebinning(rebinner)
    next_step.wrapped_step = following_sequential_step
    return next_step

def build_shear_catalog_noise(procedure_config,following_sequential_step):
    """
    Adds noise to the shear catalog.
    """

    start_seed = procedure_config.get_noise_seed()
    rs_correction = procedure_config.noise_rs_correction()
    diff_seed = procedure_config.diff_noise_seed_per_tomo_bin()
    out = NoiseAdditionStep(CatalogShapeNoiseAdder(start_seed,
                                                   diff_seed_elem = diff_seed,
                                                   rs_correction =rs_correction,
                                                   inplace = False))
    out.wrapped_step = following_sequential_step
    return out

def _build_procedure_from_shear_cat(begin_object, procedure_config,
                                    storage_collection, objects_to_save,
                                    num_tomo_bins =-1,ppz_noise=None):
    """
    Tries to build the procedure which can go as far back as converting the 
    shear catalog into the shear map.

    This function assumes that if the shear catalog is used noise has already 
    been added.
    """

    if num_tomo_bins == -1:
        tomo = False
    else:
        assert num_tomo_bins>=0
        tomo = True

    if tomo:
        extra_descr = Descriptor.tomo
    else:
        extra_descr = Descriptor.none

    # first we see if we need to build peak counting
    if _check_build_peak_counting(save_config=objects_to_save,tomo=tomo):
        feature_step = build_peak_counting(begin_object, procedure_config,
                                           storage_collection,
                                           objects_to_save,
                                           tomo = tomo,
                                           num_tomo_bins = num_tomo_bins)
        if begin_object == (extra_descr,"peak_loc"):
            return feature_step
    else:
        feature_step = None

    if tomo and objects_to_save.feature_products.tomo_power_spectrum:
        other_feat_step = build_power_spectrum(begin_object, procedure_config,
                                               storage_collection,
                                               objects_to_save,
                                               tomo = tomo)
        if feature_step is None:
            feature_step = other_feat_step
        else:
            temp = CompositeProcedureStep()
            temp.add(feature_step)
            temp.add(other_feat_step)
            feature_step = temp

    if begin_object == (Descriptor.smoothed_noisy | extra_descr, "conv_map"):
        return feature_step

    # create the step to generate the smoothed noisy convergence map
    recent_step = build_smooth_noisy_convergence(begin_object, procedure_config,
                                                 storage_collection,
                                                 objects_to_save,
                                                 feature_step, tomo = tomo)

    if begin_object == (Descriptor.noisy | extra_descr, "shear_map"):
        return recent_step

    # create the step to generate the noisy shear map
    recent_step = build_shear_conversion(begin_object, procedure_config,
                                         storage_collection, objects_to_save,
                                         following_sequential_step =recent_step,
                                         tomo = tomo)

    if recent_step is None:
        return None
    
    if begin_object == (Descriptor.none,"shear_cat"):
        if (ppz_noise is not None) and (not tomo):
            raise ValueError("Cannot run pseudo photoz process with a "
                             "non-tomographic dataset.")
        if tomo and (ppz_noise is not None 
                     or procedure_config.has_non_ppz_rebinning()):
            recent_step = build_rebinning(begin_object, procedure_config,
                                          storage_collection, objects_to_save,
                                          following_sequential_step=recent_step,
                                          num_tomo_bins = num_tomo_bins,
                                          ppz_noise = ppz_noise)
        return recent_step
    else:
        raise ValueError("Not applicable beginning step")

def _build_procedure_helper(begin_object, procedure_config,
                            storage_collection, objects_to_save,
                            num_tomo_bins,ppz_noise=None):
    # first we will try to build the non-tomographic steps
    non_tomo_proc = _build_procedure_from_shear_cat(begin_object,
                                                    procedure_config,
                                                    storage_collection,
                                                    objects_to_save,
                                                    num_tomo_bins =-1,
                                                    ppz_noise=ppz_noise)

    # next we will try to build the tomographic steps
    tomo_proc = _build_procedure_from_shear_cat(begin_object,
                                                procedure_config,
                                                storage_collection,
                                                objects_to_save,
                                                num_tomo_bins = num_tomo_bins,
                                                ppz_noise = ppz_noise)

    if tomo_proc is not None:
        # need to wrap the steps for rebinning if necessary
        if non_tomo_proc is not None:
            if begin_object != (Descriptor.none, "shear_cat"):
                raise NotImplementedError("Currently can only start tomo and "
                                          "non-tomo together from the "
                                          "shear_cat ")
            next_step = CompositeProcedureStep()
            next_step.add(non_tomo_proc)
            next_step.add(tomo_proc)
        else:
            next_step = tomo_proc
    elif non_tomo_proc is not None:
        next_step = non_tomo_proc
    else:
        raise ValueError("Unclear instructions. Instruction specify neither "
                         "computing tomographic or non-tomographic quantities")

    # for now assume we want objects with noise
    if begin_object == (Descriptor.none, "shear_cat"):
        return build_shear_catalog_noise(procedure_config, next_step)
    else:
        return next_step

def build_procedure(begin_object, procedure_config, save_config,
                    storage_collection, num_tomo_bins, ppz_noise=None):

    begin_object = (_convert_descriptor(begin_object[0]), begin_object[1])

    objects_to_save = copy.deepcopy(save_config)
    step = _build_procedure_helper(begin_object, procedure_config,
                                   storage_collection, objects_to_save,
                                   num_tomo_bins,ppz_noise)

    remaining = filter(lambda key : objects_to_save[key],
                       objects_to_save.keys())
    if len(remaining)!=0:
        
        raise ValueError(("Not all specified analysis products are scheduled "
                          "for saving. The remaining products to be saved are: "
                          "\n{:s}").format(str(remaining)))
    return step
