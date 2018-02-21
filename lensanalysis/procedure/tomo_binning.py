import numpy as np
import astropy.table as tbl

from .procedure import IntermediateProcedureStep
from ..misc.log import logprocedure


def _lazy_rebinner(orig_tomo_bins,new_intervals,colname):
    """
    We are literally just relying on the built in rebinner. We can probably be 
    smarter about this.
    
    Assumes that the new_intervals are monotomic
    """

    raise RuntimeError("I'm not convinced this works!")
    num_intervals = len(new_intervals)
    rebinned_input = [[] for i in range(num_intervals)]
    for orig_tomo_bin in orig_tomo_bins:
        values = orig_tomo_bin[colname]
        min_val = np.amin(values)
        max_val = np.amax(values)

        first_interval_index = None
        test_intervals = []
        for j,interval in enumerate(new_intervals):
            if min_val < interval[1]:
                if max_val >= interval[0]:
                    # if this is not satisfied, then the value doesn't work at
                    # all.
                    first_interval_index = j
                break

        if first_interval_index is None:
            # the current old tomo bin will not contribute to any of the new
            # tomo
            # bins
            break

        for j,interval in enumerate(new_intervals[start_index:]):
            if max_val < interval[1]:
                last_interval_index = j
                break
        else:
            last_interval_index = len(new_intervals)-1

        if first_interval_index == last_interval_index:
            rebinned_input[first_interval_index].append(orig_tomo_bin)
        else:
            start = first_interval_index
            stop = last_interval_index+1
            temp = orig_tomo_bin.rebin(new_intervals[start:stop],colname)
            for j,bin_contribution in enumerate(temp):
                rebinned_input[j+start].append(bin_contribution)

    # at this point, rebinned_input[i] contains a list of all contributions of
    # all of the old bins to the new tomographic bin i
    out = [tbl.vstack(elem) for elem in rebinned_input]

    return out

class StaticRebinner(object):
    """
    This object handles Rebinning if it does not change between realizations.

    Parameters
    ----------
    new_bin_source_loc : iterable of iterables
        Provides the location for what needs to be rebinned. More specifically, 
        the length of the iterable (lets call it the outer iterable) is the 
        number of bins we produce for this operation. Each element of the 
        outer iterable is also an iterable, (let's call it the location 
        iterable) that detail the origins of the data for a different output 
        tomographic bin. Elements of a given location iterable are themselves 
        3-element tuples (or any sequence of 3 elements). Each tuple describes 
        which elements of a particular output bin are part of the given output 
        bin. Specifically, the first entry gives the index of one of the input 
        bins, the second element gives the locations of the input bin to be 
        included in the output bin, and the third element gives the number of 
        elements that this input bin contributes to the output bin. The second 
        element can be any argument that an Astropy Column object will except 
        as indexes.
    column_data : sequence of (colname,dtype)
        Sequence of tuples listing the column name and np.dtype for each of the 
        columns to be included in the catalogs of each new bins.
    new_bin_cat_lengths : sequence of positive integers
        The lengths of each of the new bins. This must be the same length as 
        new_bin_source_loc
    share_position_component : bool, optional
        Whether or not the component of the catalog should be shared for 
        different realizations (The underlying position and redshift data would 
        be pointed to by catalogs of different realizations. While this would 
        save time and memory, any modifications to this data would affect 
        future realizations)
    """
    def __init__(self,new_bin_source_loc,column_data, new_bin_cat_lengths,
                 share_position_component=True):
        pass

    def rebin(self, data_object, map_id):
        pass

    @classmethod
    def from_old_bin_source_loc(cls, old_bin_source_loc, colnames=None,
                                share_position_component):
        pass

class DynamicRebinner(object):
    """
    This object handles rebinning that can change between realization.
    """

    def _determine_bins(self,data_object,map_id):
        """
        Determines the intervals for use while binning. Optionally modifies the 
        data_object.

        Returns
        -------
        interval: list of tuples
            List of tuples of new intervals
        mod_data_object: list of lenstools.catalog.Catalog objects
            Collection of data_objects to perform rebinning on. These may or 
            may not be modified from the input data_object
        colname : str
            The name of the column in the modified data objects used for 
            rebinning.
        """
        raise NotImplementedError()

    def rebin(self, data_object,map_id):
        interval, mod_data_object, colname = self._determine_bins(data_object,
                                                                  map_id)
        return _lazy_rebinner(mod_data_object,interval,colname)

class PhotozNoiseAddition(object):
    """
    This object is responsible for simulating photometric Errors.
    """

    def __call__(self,zspec,map_id,bin_num):
        raise NotImplementedError()

class ConstantBias(PhotozNoiseAddition):
    def __init__(self,bias):
        if bias < 0:
            raise RuntimeError("Need to come up with rules for shifts that "
                               "give negative photoz")
        self.bias = bias

    def __call__(self,zspec,map_id,bin_num):
        return self.bias*(1.+zspec)


class PseudoPhotozRebin(DynamicRebinner):
    """
    This file imitates photometric redshift errors by applying photo-z errors
    to the input catalog and then rebinning the catalog by the photo-z values.

    This method was employed for the LSST Tomography paper.

    For now this just uses a constant bias.

    Parameters
    ----------
    noise_function : PhotozNoiseAddition
        This gets called to add noise.
    bin_intervals : list of tuples, optional
        List of tuples of new intervals to rebin the catalogs to. If not 
        specified, then we will automatically use the same redshift intervals
        as were used in the input catalog.
    min_bin_zero : bool, optional
        If bin_intervals was not specified. This asks if the lowest bin 
        interval should end at z = 0 because this may not be the case. If False,
        any galaxy with photoz between z=0 and the lowest bin are discarded. 
        Default is False.
    update_in_place : bool, optional
        Whether or not the Shear Catalogs should have their redshifts updated 
        in place. Presently, this is required. The actual binning of the 
        different galaxies is unaffected.
    save_file :
        This is a placeholder - I'm not sure we will actually do anything with 
        it.
    """

    def __init__(self, noise_function,
                 bin_intervals=None,
                 min_bin_zero = False,
                 update_in_place = True,
                 save_file=None):
        self.noise_function = noise_function
        self.bin_intervals = bin_interval
        self.min_bin_zero = min_bin_zero
        assert update_in_place
        self.update_in_place = update_in_place
        assert save_file is None
        self.save_file = save_file

    def _determine_bins(self,data_object,map_id):
        func = self.noise_function

        if self.bin_interval is None:
            bin_interval = []
            # need to check whether intervals are inclusive or exclusive
            for i,catalog in enumerate(data_object):
                pass

        for i,catalog in enumerate(data_object):
            photoz = func(catalog["z"],map_id,bin_num)
            if in_place:
                catalog["z"] = photoz
            else:
                raise NotImplementedError()
        return data_object


class ShearCatRebinning(IntermediateProcedureStep):
    """
    This step rebins the shear catalogs - it usually increases or decreases 
    (but does not completely consolidate) the number of tomographic bins.

    In principle, this should also support discarding unecessary tomographic 
    information.

    For now, this is implemented with the assumption that individual bins will 
    be subdivided. In the future, if more complex rebinning is required, we can 
    refactor.
    """

    def __init__(self):
        pass
