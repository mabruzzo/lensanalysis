from operator import add

import numpy as np
import astropy.table as tbl

from lenstools.utils.algorithms import step

from .procedure import IntermediateProcedureStep
from ..misc.log import logprocedure


def _modified_step(x,intervals,vals):
    """
    Modified version of the step function defines in LensTools.

    The step function in LensTools inconsistently treats values that fall on 
    the edges of the specified intervals. In this function, we assume that the 
    intervals are monotonic. An interval is defined as a tuple: (lower,upper). 
    The range of values included in an interval are defined as [lower,upper).
    """

    # first has been defined to map np.sign(x-bound)[np.sign(x-bound)==0]=1
    # second has been defined to map np.sign(x-bound)[np.sign(bound-x)==0]=-1
    # we can probably come up with a more clever way to do this.
    first = lambda x,bound : np.where(np.sign(x-bound)>=0,1,-1)
    second = lambda x,bound : np.where(np.sign(bound-x)>0,1,-1)
    return reduce(add,(0.5*vals[n]*(first(x,i[0])+second(x,i[1]))
                       for n,i in enumerate(intervals)))

def _modified_rebin(catalog,intervals,field="z"):
    """
    This is just a modified version of the bound rebin method of 
    lenstools.catalog.Catalog. We use the exact same implementation except we 
    update the change which step function is called.
    """
    #print catalog
    #print intervals
    catalog_columns = catalog.colnames

    #Group by column interval
    catalog["group"] = _modified_step(catalog[field],intervals,
                                      np.array(range(1,len(intervals)+1)))
    catalog_rebinned = list()
    for n,i in enumerate(intervals):
        catalog_rebinned.append(catalog[catalog["group"]==n+1][catalog_columns])
    #Return re-binned catalog to user
    return catalog_rebinned

def _lazy_rebinner(orig_tomo_bins,new_intervals,colname):
    """
    We are literally just relying on the built in rebinner. We can probably be 
    smarter about this.
    
    Assumes that the new_intervals are monotomic
    """

    #raise RuntimeError("I'm not convinced this works!")
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

        for j,interval in enumerate(new_intervals[first_interval_index:]):
            if max_val < interval[1]:
                last_interval_index = j + first_interval_index
                break
        else:
            last_interval_index = len(new_intervals)

        #print "first_interval_index = {:d}".format(first_interval_index)
        #print "last_interval_index = {:d}".format(last_interval_index)
        if first_interval_index == last_interval_index:
            rebinned_input[first_interval_index].append(orig_tomo_bin)
        else:
            start = first_interval_index
            stop = last_interval_index+1
            #temp = orig_tomo_bin.rebin(new_intervals[start:stop],colname)
            temp = _modified_rebin(orig_tomo_bin,new_intervals[start:stop],
                                   colname)
            for j,bin_contribution in enumerate(temp):
                rebinned_input[j+start].append(bin_contribution)

    # at this point, rebinned_input[i] contains a list of all contributions of
    # all of the old bins to the new tomographic bin i
    #print rebinned_input
    out = []
    for i,elem in enumerate(rebinned_input):
        if len(elem) > 0:
            out.append(tbl.vstack(elem))
        else:
            raise ValueError("Bin {:d} has no entries.".format(i+1))

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
                                share_position_component=None):
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
        self.bias = bias

        self.can_be_neg = False
        if bias<0:
            self.can_be_neg = True

    def __call__(self,zspec,map_id,bin_num):
        return self.bias*(1.+zspec)

class InvertedConstantBias(PhotozNoiseAddition):
    """
    Assumes that noise photoz noise is modelled as
    bph(zs) = bias * (1+zs)

    This function assumes that we are given the photometric redshift (with 
    only a constant bias) and we want to get the spectroscopic redshift back.
    
    zp = bph + zs -> zp = bias*(1+zs) + zs -> zp = bias+(1+bias)*zs
    zs*(1+bias) = zp - bias
    zs = (zp-bias)/(1+bias)
    """
    def __init__(self,bias):
        self.bias = bias

        self.can_be_neg = False
        if bias>0:
            self.can_be_neg = True
    def __call__(self,zspec,map_id,bin_num):
        return (zspec-bias)/(1. + bias)


class PseudoPhotozRebinner(DynamicRebinner):
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
        as were used in the input catalog. The intervals are assumed to be 
        monotonically increasing. Within a tuple representing an interval, 
        the bin is assumed to be inclusive at the minimum value and exclusive 
        at the maximum value.
    colname : str,optional
        The name of the photoz column. By default it is 'z'.
    min_bin_value : float, optional
        If bin_intervals was not specified. This asks for the exclusive value 
        that the highest bin interval should end at. If this is None, then the 
        value is just set to the minimum value initially in the lowest bin.
    max_bin_value : float, optional
        If bin_intervals was not specified. This asks for the exclusive value 
        that the highest bin interval should end at. If it is not provided and 
        bin interval was not specified, the exclusive limit is set to the next 
        floating point value after the floating point value already in the bin.
    contact_intervals : bool, optional
        If bin_intervals was not specified. This asks if the upper limit of one 
        interval is always equal to the lower limit of the next interval. If 
        this is not the case, then there can be gaps between intervals. Default 
        is True. 
    update_in_place : bool, optional
        Whether or not the Shear Catalogs should have their redshifts updated 
        in place. Presently, this is required. The actual binning of the 
        different galaxies is unaffected.
    save_file :
        This is a placeholder - I'm not sure we will actually do anything with 
        it.
    """

    def __init__(self, noise_function, bin_intervals=None, colname = 'z',
                 min_bin_value = False, max_bin_value = None,
                 contact_intervals = True, update_in_place = True,
                 save_file=None):
        self.noise_function = noise_function
        self.bin_intervals = bin_intervals
        self.colname = colname

        # these next 3 variables are only important if we are dynamically
        # identifying the intervals
        self.min_bin_value = min_bin_value
        self.max_bin_value = max_bin_value
        self.contact_intervals = contact_intervals

        assert update_in_place
        self.update_in_place = update_in_place
        assert save_file is None
        self.save_file = save_file

    def _determine_bins(self,data_object,map_id):
        func = self.noise_function

        if self.bin_intervals is None:
            bin_intervals = []

            n = len(data_object)
            for i,catalog in enumerate(data_object):
                if i == 0 and self.min_bin_value is not None:
                    min_val = self.min_bin_value
                elif self.contact_intervals and i>0:
                    min_val = max_val
                else:
                    min_val = np.amin(catalog[self.colname])

                if self.contact_intervals and i<(n-1):
                    max_val = np.amin(data_object[i+1][self.colname])
                elif i == (n-1) and (self.max_bin_value is not None):
                    max_val = self.max_bin_value
                else:
                    max_val = np.nextafter(np.amax(catalog[self.colname]),
                                           np.inf)
                bin_intervals.append((min_val,max_val))
            self.bin_intervals = bin_intervals
            print self.bin_intervals

        in_place = self.update_in_place
        for i,catalog in enumerate(data_object):
            photoz = func(catalog["z"],map_id,bin_num = i+1)
            if in_place:
                catalog["z"] += photoz
                if func.can_be_neg:
                    w = (catalog["z"] <0.0)
                    catalog["z"][w] = 0.0
            else:
                raise NotImplementedError()
        return self.bin_intervals, data_object, self.colname


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

    def __init__(self,rebinner):
        self.rebinner = rebinner

    def intermediate_operation(self,data_object,packet):
        if self.rebinner is None:
            return data_object
        else:
            return self.rebinner.rebin(data_object,packet.data_id)
