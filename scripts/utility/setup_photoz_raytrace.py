"""
This script is designed to be used to:
1) Generate a new subdirectory in 
    /work/05274/tg845732/stampede2/simData/LSST100Fid/Home/photoz
   for the new cosmology
2) Write the catalog.ini file in that subdirectory for the photoz raytracing.
3) Generate the updated position_bins{:d}.fits files
4) generate the Ray.sh file & launch it/ just execute the raytracing.
"""

import argparse
import ConfigParser
import os, os.path

from lenstools.pipeline.deploy import ParsedHandler

def write_catalog_configuration(template_file,outfile,photoz_name):
    config = ConfigParser.SafeConfigParser()
    config.read(template_file)
    config.set("CatalogSettings","directory_name",photoz_name)
    with open(outfile, 'w') as config_file:
        config.write(config_file)


def write_pos_file_builder_slurm_file():
    pass


def setup(photoz_name):

    # should probably check that the photoz_name is actually contained by the
    # configuration file
    
    root_dir = '/work/05274/tg845732/stampede2/simData/LSST100Fid/Home/photoz'
    template_file = None
    new_dir = os.path.join(root_dir,photoz_name)

    if os.path.exists(new_dir):
        if not os.path.isdir(new_dir):
            raise RuntimeError(("{:s} already exists and it is not a "
                                "directory").format(new_dir))
    else:
        os.path.mkdir(new_dir)

    config_fname = os.path.join(new_dir,"catalog.ini")

    if os.path.exists(config_fname):
        if not os.path.isfile(config_fname):
            raise RuntimeError(("{:s} already exists and it is not a "
                                "file").format(config_fname))
    else:
        write_catalog_configuration(template_file,config_fname,photoz_name)

    
    # need to create build_pos_files.sh
    # need to create ray.sh


parser = argparse.ArgumentParser()
parser.add_argument("-j", "--job", dest="job_options_file", action="store",
                    type=str, required = True,
                    help="job specifications file")
parser.add_argument("-s", "--system", dest="system", action="store",
                    type=str, required = True,
                    help=("configuration file that contains the cluster "
                          "specifications"))
parser.add_argument("-r", "--root_dir", dest="root_dir", action="store",
                    type=str, required = True
                    help=("root directory within which we construct the new "
                          "directory for a given photoz scheme and populate "
                          "with the required configuration files."))
parser.add_argument("-c", "--config", dest = "config", required = True,
                    help = ("Specify the photoz file that includes information "
                            "about different photoz simulation schemes."))
parser.add_argument("-t", "--template_cat", dest = "template_cat",
                    required = True,
                    help = "Specify the path to the template catalog file.")
parser.add_argument("id",nargs=1, required = True,
                    help = "The name of the photoz error scheme to use.")



if __name__ == '__main__':
    write_catalog_configuration("catalog.ini","test_cat.ini","posConstBias")
