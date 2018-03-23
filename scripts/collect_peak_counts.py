import argparse
import multiprocessing
import os.path
from string import Formatter
import sys

import numpy as np

from lensanalysis.config_parsing.cosmo_config import \
    CosmoCollectionConfigBuilder
from lensanalysis.misc.analysis_collection import default_value_UAPC

"""
This script is designed to collect all of the computed peak counts.
"""


parser = argparse.ArgumentParser()

parser.add_argument("--name",dest = "names", nargs='+', default = [],
                    action = 'store',
                    help = ("The name of cosmology to post-process. If this is "
                            "none, then it processes all cosmologies"))

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


def collector(fname, loader, start, stop, average = False, tomo = True,
              num_tomo_bins = 0):

    if tomo:
        temp = [[] for i in range(num_tomo_bins)]
    else:
        temp = []

    for i in xrange(start,stop):
        try:
            peak_count_col = loader.load(i)
        except ValueError:
            print "Error occured for realization {:d}".format(i)
            print "Was supposed to be saved in: {:s}".format(fname)
            raise 
        if tomo:
            for j,peak_counts in enumerate(peak_count_col):
                # check to make sure we get the right shape
                temp[j].append(peak_counts._counts)
        else:
            temp.append(peak_count_col[0]._counts)

    if tomo:
        out = []
        for collection in temp:
            if average:
                out.append(np.mean(np.stack(collection),axis = 0))
            else:
                out.append(np.stack(collection))
        result = np.stack(out)
    else:
        if average:
            result = np.average(np.stack(temp),axis=0)
        else:
            result = np.stack(temp)
    print result.shape
    np.save(fname,result)

def process_name(name, fname_template, cosmo_storage_col, load_config, start,
                 stop, tomo, cmd_args, num_tomo_bins):
    if check_num_fields(fname_template) == 0:
        fname = fname_template
    else:
        fname = fname_template.format(name)

    storage = cosmo_storage_col.get_analysis_product_storage(name,
                                                             load_config)
    if tomo:
        loader = storage.feature_products.tomo_peak_counts
    else:
        loader = storage.feature_products.peak_counts

    collector(fname, loader, start, stop, cmd_args.average, tomo = tomo,
              num_tomo_bins = num_tomo_bins)

class ProcessNameWrapper(object):
    def __init__(self, fname_template, cosmo_storage_col, load_config, start,
                 stop, tomo, cmd_args, num_tomo_bins):
        self.fname_template = fname_template 
        self.cosmo_storage_col = cosmo_storage_col
        self.load_config = load_config 
        self.start = start
        self.stop = stop
        self.tomo = tomo
        self.cmd_args = cmd_args
        self.num_tomo_bins = num_tomo_bins

    def __call__(self,name):
        process_name(name, self.fname_template, self.cosmo_storage_col, 
                     self.load_config, self.start, self.stop, self.tomo, 
                     self.cmd_args, self.num_tomo_bins)

if __name__ == '__main__':
    print sys.argv
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
    tomo = cmd_args.tomo
    load_config = default_value_UAPC(False)
    if tomo:
        load_config.feature_products.tomo_peak_counts = True
        num_tomo_bins = cosmo_storage_col.get_num_tomo_bin()
    else:
        load_config.feature_products.peak_counts = True
        num_tomo_bins = 0

    if cmd_args.mp_collections is None:
        num_procs = 1
    else:
        num_procs = cmd_args.mp_collections
        assert num_procs >=1 and isinstance(num_procs,int)

    f = ProcessNameWrapper(fname_template, cosmo_storage_col, load_config,
                           start, stop, tomo, cmd_args, num_tomo_bins)
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
