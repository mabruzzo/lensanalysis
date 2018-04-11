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
from lenstools.pipeline.settings import JobSettings, CatalogSettings

from add_photoz_errors import get_input_template_and_indices,get_indices

def write_catalog_configuration(template_file,outfile,photoz_name):
    config = ConfigParser.SafeConfigParser()
    config.read(template_file)
    config.set("CatalogSettings","directory_name",photoz_name)
    with open(outfile, 'w') as config_file:
        config.write(config_file)

def update_stdouterr_redirection(job_settings,script_filename):
    # update the redirections of stdout and stderr to be in the directory where
    # the script has been placed.

    # get directory where script is being saved
    out_dir = os.path.basename(os.path.abspath(script_filename))

    for attr in ['redirect_stdout','redirect_stderr']:
        initial = getattr(job_settings,attr)
        basename, fname = os.path.split(initial)
        if basename != '':
            raise ValueError(("The {:s} option must not include a directory "
                              "path").format(attr))
        setattr(job_settings,attr,os.path.join(out_dir,fname))

def write_pos_file_builder_slurm_script(path,job_file,job_handler,
                                        input_template,start,stop,
                                        photoz_config_fname,
                                        photoz_name):
    script_filename = path

    job_settings = JobSettings.read(job_file, "PosFileBuilding")
    update_stdouterr_redirection(job_settings,script_filename)

    num_proc = job_settings.cores_per_simulation
    arguments = ["-i {:s}".format(os.path.abspath(input_template)),
                 "-o position_cat{:d}.fits",
                 "-b {:d}".format(start),
                 "-s {:d}".format(stop),
                 "-c {:s}".format(os.path.abspath(photoz_config_fname))]
    if num_proc > 1:
        arguments.append("--mp {:d}".format(num_proc))
    arguments.append("{:s}".format(photoz_name))

    temp = ["python", job_settings.path_to_executable]
    executable = ' '.join(temp+arguments)
    working_dir = os.path.dirname(os.path.abspath(path))

    job_settings.num_cores = job_settings.cores_per_simulation
    with open(script_filename,'w') as scriptfile:
        scriptfile.write(job_handler.writePreamble(job_settings))
        scriptfile.write("cd {:s}\n".format(working_dir))
        scriptfile.write('{:s} > building.log'.format(executable))
        scriptfile.write('\n')


def write_raytrace_script(path,job_file,job_handler,ic_id):
    # need to do this in one simple chunk
    # use the lensplanes saved under the name of the fiducial cosmology

    """ 
    This is basically a modified version of the writeRaySubmission bound 
    method defined in the SimulationBatch class of LensTools
    """

    script_filename = path
    job_settings = JobSettings.read(job_file, "RayTracing")

    update_stdouterr_redirection(job_settings,script_filename)

    job_settings.num_cores = job_settings.cores_per_simulation

    parts = ic_id.split("|")

    if len(parts) == 2:
        try:
            config_fname = os.path.join(os.path.dirname(path),"catalog.ini")

            raytracing_settings = CatalogSettings.read(config_fname)
            if (raytracing_settings.lens_catalog_realizations
                % job_settings.cores_per_simulation):
		raise ValueError("The number of map realizations must be a "
                                 "multiple of the number of cores per "
                                 "simulation!")
        except AssertionError:
            raise NotImplementedError("Not currently able to handle Planes of "
                                      "source galaxies.")
    elif len(parts) == 1:
        raise NotImplementedError("Not curently able handle "
                                  "TelescopicMapSettings.")
    else:
        raise ValueError(("There are too many '|'' in your id: "
                          "{0}".format(realizations_in_chunk[e])))
    
    executable = (job_settings.path_to_executable + " " +
                  '-e {0} -c {1} "{2}" '.format(environment_file, config_fname,
                                                ic_id))
    with self.syshandler.open(script_filename,"w") as scriptfile:
        scriptfile.write(job_handler.writePreamble(job_settings))
        scriptfile.write("cd {:s}\n".format(working_dir))
        scriptfile.write(job_handler.writeExecution([executable],
                                                    job_settings.num_cores,
                                                    job_settings))

def setup(photoz_name,system_file, job_file, catalog_template, input_template,
          start, stop, photoz_config_fname,ic_id):

    # should probably check that the photoz_name is actually contained by the
    # configuration file

    job_handler = ParsedHandler.read(system_file)
    
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
        write_catalog_configuration(catalog_template,config_fname,photoz_name)

    # create build_pos_files.sh
    pos_file_script = os.path.abspath(os.path.join(new_dir,
                                                   "build_pos_files.sh"))
    if os.path.exists(pos_file_script):
        if not os.path.isfile(pos_file_fname):
            raise RuntimeError(("{:s} already exists and it is not a "
                                "file").format(pos_file_fname))
    else:
        write_pos_file_builder_slurm_script(pos_file_script, job_file,
                                            job_handler, input_template,
                                            start, stop, photoz_config_fname,
                                            photoz_name)


    # create ray.sh
    ray_file_script = os.path.abspath(os.path.join(new_dir,"ray.sh"))
    if os.path.exists(ray_file_script):
        if not os.path.isfile(ray_file_fname):
            raise RuntimeError(("{:s} already exists and it is not a "
                                "file").format(ray_file_fname))
    else:
        # need to create ray.sh
        write_raytrace_script(ray_file_script,job_file,job_handler,ic_id)

def get_config_file(cmd_args,attr_name,option_name):
    fname = getattr(cmd_args,attr_name)
    if os.path.isfile(fname):
        return fname
    raise ValueError(("The file name specified for {:s}, {:s}, does not point "
                      "to an existing file.").format(option_name,fname))

def driver(cmd_args):
    start,stop = get_indices(cmd_args)
    input_template = get_input_template_and_indices(cmd_args)[0]
    job_file = get_config_file(cmd_args,"job_options_file","--job")
    system_file = get_config_file(cmd_args,"system","--system")
    photoz_config_fname = get_config_file(cmd_args,"config","--config")
    catalog_template = get_config_file(cmd_args,"template_cat",
                                       "--template_cat")

    ic_id = cmd_args.model_id
    photoz_name = cmd_args.id
    setup(photoz_name,system_file, job_file, catalog_template,
          input_template, start, stop, photoz_config_fname, ic_id)



parser = argparse.ArgumentParser()
parser.add_argument("-j", "--job", dest="job_options_file", action="store",
                    type=str, required = True,
                    help="job specifications file")
parser.add_argument("-s", "--system", dest="system", action="store",
                    type=str, required = True,
                    help=("configuration file that contains the cluster "
                          "specifications"))
parser.add_argument("-r", "--root_dir", dest="root_dir", action="store",
                    type=str, required = True,
                    help=("root directory within which we construct the new "
                          "directory for a given photoz scheme and populate "
                          "with the required configuration files."))
parser.add_argument("-c", "--config", dest = "config", required = True,
                    help = ("Specify the photoz file that includes information "
                            "about different photoz simulation schemes."))
parser.add_argument("-t", "--template_cat", dest = "template_cat",
                    required = True,
                    help = "Specify the path to the template catalog file.")
parser.add_argument("id",nargs=1,
                    help = "The name of the photoz error scheme to use.")
parser.add_argument("-i","--input_template", dest = "input_template",
                    action = "store", default = None, required = True,
                    help = ("The template with which the input file(s) should "
                            "be saved. This template should include 1 "
                            "formatter option corresponding to where the "
                            "index. Example: .../.../position_cat{:d}.fits"))
parser.add_argument("-b", "--begin", dest = "begin", action = "store",
                    type = int, required = True,
                    help = ("Specify the minimum index to use for formatting "
                            "the templates"))
parser.add_argument("--stop", dest = "stop", action = "store",
                    type = int, required = True,
                    help = ("Specify the minimum index to use for formatting "
                            "the templates ( (--stop - 1) gives the final "
                            "index used to format a template)."))
parser.add_argument("-m","--model_id", dest = "model_id", required = True,
                    help = ('Supplies a file which contains the ic of '
                            'the planes to use for raytracing, in the form '
                            '"cosmo_id|geometry_id". '
                            'E.g: Om0.260_Ol0.740_w-1.000_si0.800|512b260'))


if __name__ == '__main__':
    cmd_args = parser.parse_args()
    driver(cmd_args)
    """
    job_handler = ParsedHandler.read("../../lensanalysis/sample_config/"
                                     "utility_config/stampede2.ini")
    write_pos_file_builder_slurm_script('build_pos_files.sh',
                                        ('../../lensanalysis/sample_config/'
                                         'utility_config/setup_photoz_job.ini'),
                                        job_handler,
                                        '../../jobs/position_cat{:d}.fits',
                                        1,6,
                                        ('../../lensanalysis/sample_config/'
                                         'photoz_config.ini'),
                                        'constPosBias')

    write_raytrace_script(os.path.abspath("ray.sh"),
                          ('../../lensanalysis/sample_config/'
                           'utility_config/setup_photoz_job.ini'),
                          job_handler,
                          'Om0.260_Ol0.740_w-1.000_si0.800|512b260')
    """
