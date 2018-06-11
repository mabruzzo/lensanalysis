import argparse
import multiprocessing
import os.path
from string import Formatter
import sys

import numpy as np

from lensanalysis.config_parsing.cosmo_config import \
    CosmoCollectionConfigBuilder
from lensanalysis.misc.analysis_collection import default_value_UAPC
from lensanalysis.misc.enum_definitions import Descriptor

"""
This script is designed to collect all of the computed peak counts.
"""


parser = argparse.ArgumentParser()

parser.add_argument("--name",dest = "names", nargs='+', default = [],
                    action = 'store',
                    help = ("The name of cosmology to post-process. If this is "
                            "none, then it processes all cosmologies"))
parser.add_argument("-f","--feature",dest="feature", action = "store",
                    required = True,
                    help = ("Indicates the type of feature we want to collect. "
                            "Currenty allowed values include: 'peak_counts' "
                            "and 'power_spectrum'"))
parser.add_argument("--fiducial",dest = "fiducial", action = "store_true",
                    default = False,
                    help = ("Specify that the analysis is being performed on "
                            "the fiducial cosmology (rather than a sampled "
                            "cosmology."))
parser.add_argument("-c", "--config", dest = "config", required = True,
                    help = ("Specify the configuration file that details how "
                            "all files are (or should be) saved."))

parser.add_argument("--min",dest = "min_realization", action = "store",
                    type = int, required = True,
                    help = ("Specify the minimum realization to process."
                            "Default is the maximum allowed realization to "
                            "be processed."))

parser.add_argument("--max",dest = "max_realization", action = "store",
                    type =int, required = True,
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed."))
parser.add_argument("--template", dest = "template", action = "store",
                    default = None, required = True,
                    help = ("The template with which the output file(s) should "
                            "be saved."))
parser.add_argument("--tomo", dest = "tomo", action = "store_true",
                    default = False,
                    help = ("Indicates if the peak counts are tomographic."))
parser.add_argument("--average", dest = "average", action = "store_true",
                    default = False,
                    help = ("Indicates if the peak counts should be averaged."))
parser.add_argument("--mp_collections", dest = "mp_collections", 
                    action = "store", default = None, type = int,
                    help = ("Indicates if we should use multiprocessing for "
                            "collecting different peak_count collections.\n"
                            "This expects an integer number of processes "
                            "greater than 1."))


def check_num_fields(template):
    iterable = Formatter().parse(template)

    num = 0
    for _, field_name, format_spec, _ in iterable:
        if field_name is not None:
            num +=1
    return num

def check_template(template,num_names):
    if num_names == 1:
        assert check_num_fields(cmd_args.template) in [0,1]
    else:
        assert check_num_fields(cmd_args.template) in [1]
    dir_name = os.path.dirname(template)
    assert os.path.isdir(dir_name)



def _consolidate_helper(start, stop, loader, feature_retriever):
    temp = []
    assert stop - start > 1
    for i in xrange(start,stop):
        try:
            feature_col = loader.load(i)
        except ValueError:
            print "Error occured for realization {:d}".format(i)
            raise
        temp.append(feature_retriever(feature_col))
    out = np.array(temp)

    if len(out.shape) == 2:
        # this is non-tomographic
        # we need to add a third axis
        return np.array([out])
    elif len(out.shape) == 3:
        # if we are using single file storage, then the array we have
        # constructed has different tomographic bins along axis 1 and
        # different realizations along axis 0. These need to be swapped
        # so that we return the expected result
        return np.swapaxes(out,0,1)
    else:
        raise RuntimeError("All results should be 2D or 3D")
    
class Consolidator(object):
    """
    This object is designed to consolidate feature sets into a stacked array of 
    3 Dimensions

    Different entries along axis 0 correspond to data from different 
    tomographic bins (or bin combinations).
    Different entries along axis 1 corresponds to different Realizations.
    Different entires along axis 2 correspond to components to components of a 
    data vector along within a given tomographic bin (or bin combination).


    I am unclear on why I had originally made this. The way the loaders work,
    seems to alleviate the need for creating this and subclassing it.
    """

    def __init__(self, attr_name):
        self.attr_name = attr_name

    def consolidate(self,loader,start,stop):
        raise NotImplementedError()

class SingleFileConsolidator(Consolidator):
    """
    This handles feature object collections saved to a single file.

    Example: Tomographic Power Spectrum
    """

    def consolidate(self,loader,start,stop):
        attr_name = self.attr_name
        #print loader.load(3)
        #func = lambda feature_col : getattr(feature_col[0], attr_name)
        func = lambda feature_col : [getattr(elem, attr_name) for elem in \
                                     feature_col]
        return _consolidate_helper(start, stop, loader, func)
        

class FileGroupConsolidator(Consolidator):
    """
    This handles feature object collections saved to groups of files.

    Example: Peak Count histogram
    """

    def consolidate(self,loader,start,stop):
        attr_name = self.attr_name
        func = lambda feature_col : [getattr(elem, attr_name) for elem in \
                                     feature_col]
        return _consolidate_helper(start, stop, loader, func)

# Entries of the dictionary are 3 elements
# 1st - If the feature allows non-tomographic versions
# 2nd - If the feature allows tomographic versions
# 3rd - An initialized Consolidator object
valid_features = {"power_spectrum" : (None,'tomo_power_spectrum',
                                      SingleFileConsolidator('_power')),
                  "peak_counts" : ('peak_counts','tomo_peak_counts',
                                   SingleFileConsolidator('_counts'))}

def get_feature_name(cmd_args):
    feature_name = cmd_args.feature
    assert feature_name in valid_features

    if cmd_args.tomo:
        if valid_features[feature_name][1] is None:
            raise ValueError(("functionallity for {:s} is not defined for "
                              "tomography").format(feature_name))
    else:
        if not valid_features[feature_name][0] is None:
            raise ValueError(("functionallity for {:s} is not defined for "
                              "non-tomography").format(feature_name))
    return feature_name,valid_features[feature_name][2]
                    
def setup_load_config(feature_name, cosmo_storage_col):
    tomo = cmd_args.tomo
    load_config = default_value_UAPC(False)
    if tomo:
        load_config[(Descriptor.tomo,feature_name)] = True
        num_tomo_bins = cosmo_storage_col.get_num_tomo_bin()
    else:
        load_config[(Descriptor.none,feature_name)] = True
        num_tomo_bins = 0
    return tomo, load_config, num_tomo_bins

def collector(fname, loader, start, stop, consolidator, average = False,
              tomo = True):
    try:
        temp = consolidator.consolidate(loader, start, stop)
    except ValueError:
        print "Was supposed to be saved in: {:s}".format(fname)
        raise

    if average:
        out = np.mean(temp,axis=1)
    else:
        out = temp
    if not tomo:
        out = out[0,...]
    print out.shape
    np.save(fname,out)

def process_name(name, fname_template, cosmo_storage_col, load_config, start,
                 stop, tomo, cmd_args, num_tomo_bins, consolidator,
                 feature_name):
    if check_num_fields(fname_template) == 0:
        fname = fname_template
    else:
        fname = fname_template.format(name)

    storage = cosmo_storage_col.get_analysis_product_storage(name,
                                                             load_config)
    if tomo:
        attr_name = valid_features[feature_name][1]
    else:
        attr_name = valid_features[feature_name][0]
    loader = getattr(storage.feature_products, attr_name)

    collector(fname, loader, start, stop, consolidator, cmd_args.average,
              tomo = tomo)

class ProcessNameWrapper(object):
    def __init__(self, fname_template, cosmo_storage_col, load_config, start,
                 stop, tomo, cmd_args, num_tomo_bins, consolidator,
                 feature_name):
        self.fname_template = fname_template 
        self.cosmo_storage_col = cosmo_storage_col
        self.load_config = load_config 
        self.start = start
        self.stop = stop
        self.tomo = tomo
        self.cmd_args = cmd_args
        self.num_tomo_bins = num_tomo_bins
        self.consolidator = consolidator
        self.feature_name = feature_name

    def __call__(self,name):
        process_name(name, self.fname_template, self.cosmo_storage_col, 
                     self.load_config, self.start, self.stop, self.tomo, 
                     self.cmd_args, self.num_tomo_bins, self.consolidator,
                     self.feature_name)



        
if __name__ == '__main__':
    cmd_args = parser.parse_args()
    print cmd_args.names
    num_names = len(cmd_args.names)
    names = cmd_args.names
    
    assert 1<=cmd_args.min_realization <=cmd_args.max_realization
    start = cmd_args.min_realization
    stop = cmd_args.max_realization +1
    fid = cmd_args.fiducial

    builder = CosmoCollectionConfigBuilder.from_config_file(cmd_args.config)

    if cmd_args.fiducial:
        cosmo_storage_col = builder.get_fiducial_storage_collection()
    else:
        cosmo_storage_col = builder.get_sampled_storage_collection()

    if num_names == 0:
        names = cosmo_storage_col.list_analysis_product_names()
        num_names = len(names)
        assert num_names >0
    fname_template = cmd_args.template
    check_template(fname_template,num_names)
    feature_name, consolidator = get_feature_name(cmd_args)
    tomo, load_config, num_tomo_bins = setup_load_config(feature_name,
                                                         cosmo_storage_col)

    if cmd_args.mp_collections is None:
        num_procs = 1
    else:
        num_procs = cmd_args.mp_collections
        assert num_procs >=1 and isinstance(num_procs,int)

    f = ProcessNameWrapper(fname_template, cosmo_storage_col, load_config,
                           start, stop, tomo, cmd_args, num_tomo_bins,
                           consolidator, feature_name)
    if num_procs == 1 or len(names) == 1:
        map(f,names)
    else:
        
        if len(names) <= num_procs:
            chunk_size=1
            p = multiprocessing.Pool(len(names))
        else:
            chunk_size = num_procs//num_procs
            p = multiprocessing.Pool(num_procs)
        p.map(f,names,chunk_size)
