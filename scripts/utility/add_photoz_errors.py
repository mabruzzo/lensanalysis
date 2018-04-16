import argparse
from multiprocessing import Pool
import os.path
from string import Formatter

import numpy as np

from lenstools.catalog.shear import Catalog
from lensanalysis.config_parsing.photoz_config import PhotozConfig

description = ("A program that applies photoz effects to position catalogs. It "
               "is meant to be used so that raytracing can be run on the "
               "output.")

parser = argparse.ArgumentParser(description = description)

parser.add_argument("-i", "--input_template", dest = "input_template",
                    action = "store", default = None, required = True,
                    help = ("The template with which the input file(s) should "
                            "be saved. This template should include 1 "
                            "formatter option corresponding to where the "
                            "index. Example: .../.../position_cat{:d}.fits"))
parser.add_argument("-o", "--output_template", dest = "output_template",
                    action = "store", default = None, required = True,
                    help = ("The template with which the output file(s) should "
                            "be saved. This template should include 1 "
                            "formatter option corresponding to where the "
                            "index. Example: .../.../position_cat{:d}.fits"))
parser.add_argument("-b", "--begin", dest = "begin", action = "store",
                    type = int, required = True,
                    help = ("Specify the minimum index to use for formatting "
                            "the templates"))
parser.add_argument("-s", "--stop", dest = "stop", action = "store",
                    type = int, required = True,
                    help = ("Specify the minimum index to use for formatting "
                            "the templates ( (--stop - 1) gives the final "
                            "index used to format a template)."))
parser.add_argument("-c", "--config", dest = "config", required = True,
                    help = ("Specify the photoz file that includes information "
                            "about different photoz simulation schemes."))
parser.add_argument("--mp", dest = "mp", 
                    action = "store", default = None, type = int,
                    help = ("Indicates if we should use multiprocessing for "
                            "adding photoz errors.This expects an integer "
                            "number of processes greater than 1."))
parser.add_argument("id",nargs=1,
                    help = "The name of the photoz error scheme to use.")


def get_indices(cmd_args):
    start = cmd_args.begin
    stop = cmd_args.stop

    if stop> start:
        return start,stop

    raise ValueError("The --stop argument must be larger than the --begin "
                     "argument")

def check_num_fields(template):
    """
    Copied from collect_peak_counts.py
    """
    iterable = Formatter().parse(template)

    num = 0
    for _, field_name, format_spec, _ in iterable:
        if field_name is not None:
            num +=1
    return num

def check_template(template_name,template_val):
    # check that the template has 1 field
    if check_num_fields(template_val) != 1:
        raise ValueError(('The provided value for {:s} "{:s}" must '
                          'have one field to specify an '
                          'index.').format(template_name, template_val))

    # check that the template is a fits file
    if template_val[-5:] != ".fits":
        raise ValueError(('The {:s} option must supply templates for ".fits" '
                          'files').format(template_name))

def get_input_template_and_indices(cmd_args):
    # first load the indices
    start,stop = get_indices(cmd_args)

    # next load the input template
    input_template = cmd_args.input_template

    # check the formatting of the template
    check_template("--input_template", input_template)

    # next check that all of the files exist
    for i in range(start,stop):
        path = input_template.format(i)
        if not os.path.isfile(path):
            raise ValueError("There is no file called {:s}".format(path))
    return input_template, range(start,stop)


def get_output_template(cmd_args,attr='output_template'):
    """
    repurposed from collect_peak_counts.
    """
    output_template = getattr(cmd_args,attr)

    # check the formatting of the template
    check_template("--output_template", output_template)

    dirname = os.path.dirname(output_template)
    # ensure that the output template points to files in an existing directory
    if dirname!= '' and not os.path.isdir(dirname):
        raise ValueError(('The directory in which output files formated with '
                          '--output_template will be saved, "{:s}" does not '
                          'exist').format(os.path.dirname(output_template)))
    return output_template

def get_num_processes_and_chunksize(cmd_args,num_indices):
    """
    Returns nproc, chunksize. If nproc ==1 was not specified, then it returns 
    1,None

    """
    nproc = cmd_args.mp

    if nproc is None:
        nproc = 1
        if nproc<=1:
            raise ValueError("--mp must be a positive integer greater than 1.")
        
    if nproc == 1 or num_indices ==1 :
        chunksize = None
    else:
        if num_indices <=nproc:
            nproc = num_indices
            chunksize = 1
        else:
            chunksize = num_indices // nproc
    return nproc,chunksize

def get_photoz_func(cmd_args):
    config_fname = cmd_args.config

    if not os.path.isfile(config_fname):
        raise ValueError(("The file provided for --config, {:s}, does not "
                          "exist").format(config_fname))
    config = PhotozConfig.from_fname(config_fname)

    photoz_func = config.get_photoz_noise_addition(cmd_args.id[0])
    return photoz_func

class ConfigFunc(object):
    def __init__(self,photoz_func, input_template,output_template):
        self.photoz_func = photoz_func
        self.input_template = input_template
        self.output_template = output_template

    def __call__(self,index):
        input_fname = self.input_template.format(index)
        output_fname = self.output_template.format(index)
        photoz_func = self.photoz_func

        cat = Catalog.read(input_fname)
        cat["z"] = photoz_func(cat["z"],0,index)
        if photoz_func.can_be_neg:
            w = (cat["z"] <0.0)
            cat["z"][w] = 0.0

        cat.write(output_fname)

def driver(cmd_args):
    """
    The main driver of the program.
    """

    input_template, indices = get_input_template_and_indices(cmd_args)
    output_template = get_output_template(cmd_args)
    nproc,chunksize = get_num_processes_and_chunksize(cmd_args,len(indices))

    photoz_func = get_photoz_func(cmd_args)
    func = ConfigFunc(photoz_func, input_template, output_template)

    if chunksize is None:
        map(func,indices)
    else:
        p = Pool(nproc)
        p.map(func,indices,chunksize)

    

    
if __name__ == '__main__':
    cmd_args = parser.parse_args()
    driver(cmd_args)
