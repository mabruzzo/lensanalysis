import argparse
import os.path
from string import Formatter

import numpy as np

from lensanalysis.config_parsing.cosmo_config import \
    CosmoCollectionConfigBuilder
from lensanalysis.misc.analysis_collection import default_value_UAPC

"""
This script is designed to collect all of the computed peak counts.
"""


parser = argparse.ArgumentParser()

parser.add_argument("--name", nargs='?', default = [],
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
                    type = int, required = True
                    help = ("Specify the minimum realization to process."
                            "Default is the maximum allowed realization to "
                            "be processed."))

parser.add_argument("--max",dest = "max_realization", action = "store",
                    type =int, required = True
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed."))

parser.add_argument("--template", dest = "template", action = "store",
                    default = None, required = True
                    help = ("The template with which the output file(s) should "
                            "be saved."))

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
        peak_count_col = loader.load(i)

        if tomo:
            for peak_counts, tomo_cont in zip(peak_count_cot,temp):
                tomo_cont.append(peak_counts._counts)
        else:
            temp.append(peak_count_col[0]._counts)

    if tomo:
        out = []
        for collection in temp:
            if average:
                out.append(np.mean(np.stack(temp),axis = -1))
                raise RuntimeError("make sure this works right")
            else:
                out.append(np.stack(temp))
        result = np.stack(out)
    else:
        raise RuntimeError()
    np.save(fname,result)
    raise RuntimeError("FINALLY WE NEED TO SAVE THE FILE.")


if __name__ == '__main__':
    cmd_args = parser.parse_args()
    if cmd_args.names is None:
        num_names = 0
    else:
        num_names = len(cmd_args.names)
        names = cmd_args.names
    fname_template = cmd_args.template
    check_template(fname_template,num_names)
    assert 1<=cmd_args.min_realization <=cmd_args.max_realization
    start = cmd_args.min_realization
    stop = cmd_arg.max_realization +1
    fid = cmd_args.fiducial

    builder = CosmoCollectionConfigBuilder.from_config_file(cmd_args.config)

    tomo = True
    load_config = default_value_UAPC(False)
    if tomo:
        load_config.feature_products.tomo_peak_counts = True
    else:
        load_config.feature_products.peak_counts = True

    if cmd_args.fiducial:
        cosmo_storage_col = builder.get_fiducial_storage_collection()
    else:
        cosmo_storage_col = builder.get_sampled_storage_collection()

    if num_names == 0:
        names = cosmo_storage_col.list_analysis_product_names()

    for name in names:
        if check_num_fields(template) == 0:
            fname = template
        else:
            fname = template.format(name)

        storage = cosmo_storage_col.get_analysis_product_storage(name,
                                                                 load_config)
        if tomo:
            loader = storage.feature_products.tomo_peak_counts
        else:
            loader = storage.feature_products.peak_counts

        pass
