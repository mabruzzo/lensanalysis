import subprocess
import sys
import os.path

"""
This is more for convenience. The actual script must be modified for a given 
configuration.

When setup properly, just call the script followed by the name of the 
photozError model
"""

setup_script_loc = ("/work/05274/tg845732/stampede2/lensanalysis/scripts/"
                    "utility/setup_photoz_raytrace.py")

job_file = ('/work/05274/tg845732/stampede2/raytrace_photoz_settings/'
            'setup_photoz_job.ini')
system_file = ('/work/05274/tg845732/stampede2/raytrace_photoz_settings/'
               'stampede2.ini')
root_dir = '/work/05274/tg845732/stampede2/simData/LSST100Fid/Home/photoz'
config_file = ('/work/05274/tg845732/stampede2/raytrace_photoz_settings/'
               'photoz_config.ini')
#template_cat = 'positions_bin{:d}.fits'
template_cat = ('/work/05274/tg845732/stampede2/simData/LSST100Fid/Home/'
                'catalog.ini')
input_template = ("/work/05274/tg845732/stampede2/simData/LSST100Fid/"
                  "Home/Jobs/positions_bin{:d}.fits")
#model_file = ('/work/05274/tg845732/stampede2/raytrace_photoz_settings/'
#              'collections.txt')
model_id = "'Om0.260_Ol0.740_w-1.000_si0.800|512b260'"
env_file = ("/work/05274/tg845732/stampede2/simData/LSST100Fid/Home/"
            "environment.ini")

# the indices to use with position bin templates
begin = '1'
stop = '6'

def process_output(output):
    l = output.split('\n')
    if len(l) != 18:
        print "Job submission failed. Submission output shown below:"
        print output
        raise RuntimeError()
    if l[16][:20] != "Submitted batch job ":
        print "Job submission output unexpected. Submission output shown below."
        print output
        raise RuntimeError()
    return int(l[16][20:])

if __name__ == '__main__':
    assert len(sys.argv) == 2

    photoz_name = sys.argv[1]

    # first we run the setup_photoz_raytrace.py script
    setup_command = ["python", setup_script_loc,
                     "--job",job_file,
                     "--system", system_file,
                     "--root_dir", root_dir,
                     "--config", config_file,
                     "--template_cat", template_cat,
                     "--input_template", input_template,
                     "--model_id", model_id,
                     "-e", env_file,
                     "--begin",begin, "--stop", stop,
                     photoz_name]
    print "RUNNING SETUP\n"

    code = subprocess.call(' '.join(setup_command),shell=True, 
                           stdout=sys.stdout, stderr = sys.stderr)

    # Second, check exit code of setup_photoz_raytrace.py
    # if the exit code is anything other than 0 we exit
    if code != 0:
        raise RuntimeError("setup_photoz_raytrace.py was unsuccesful")

    #raise RuntimeError()
    # Next we run add_photoz_errors.py (OPTIONALLY - may want to check that the
    # position files have not already been generated) by submitting the
    # build_pos_files.sh script

    script_loc = "{:s}/{:s}/build_pos_files.sh".format(root_dir,photoz_name)
    noise_command = ["sbatch", os.path.abspath(script_loc)]

    print "Submitting job to add photoz Errors"
    noise_output = subprocess.check_output(' '.join(noise_command),
                                           shell = True)
    noise_addition_id = process_output(noise_output)


    # Then, we need to execute the ray.sh
    script_loc = "{:s}/{:s}/ray.sh".format(root_dir,photoz_name)
    raytrace_command = ["sbatch",
                        "--dependency=afterany:{:d}".format(noise_addition_id),
                        os.path.abspath(script_loc)]
    print "Submitting Raytracing job"
    #raise RuntimeError("WE MAY WANT JOBID FOR RAYTRACING SLURM JOB")
    ray_output = subprocess.check_output(noise_command, shell = True)
    ray_id = process_output(ray_output)

    # Finally - OPTIONAL - we may want to build a score database for the
    # resulting Shear Catalogs.

