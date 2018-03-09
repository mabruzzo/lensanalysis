import subprocess
import sys
import os.path


root_dir = '/work/05274/tg845732/stampede2/simData/LSST100Fid/Home/photoz'


if __name__ == '__main__':
    assert len(sys.argv) == 2

    photoz_name = sys.argv[1]

    # first we run the setup_photoz_raytrace.py script
    setup_command = ["python", "setup_photoz_raytrace.py"]
    subprocess.call(setup_command,shell=True, stdout=sys.stdout,
                    stderr = sys.stderr)

    # Second, check exit code of setup_photoz_raytrace.py
    # if the exit code is anything other than 0 we exit



    # Next we run add_photoz_errors.py (OPTIONALLY - may want to check that the
    # position files have not already been generated) by submitting the
    # build_pos_files.sh script

    script_loc = "{:s}/{:s}/build_pos_files.sh".format(root_dir,photoz_name)
    noise_command = ["sbatch", os.path.abspath(script_loc)]

    raise RuntimeError("WE NEED TO GET THE JOBID FOR THE SUBMITTED SLURM JOB")
    subprocess.call(noise_command,shell = True, stdout=sys.stdout,
                    stderr = sys.stderr)

    noise_addition_id = None
    

    # Then, we need to execute the ray.sh
    script_loc = "{:s}/{:s}/ray.sh".format(root_dir,photoz_name)
    raytrace_command = ["sbatch",
                        "--dependency=afterany:{:d}".format(noise_addition_id),
                        os.path.abspath(script_loc)]

    #raise RuntimeError("WE MAY WANT JOBID FOR RAYTRACING SLURM JOB")
    subprocess.call(noise_command,shell = True, stdout=sys.stdout,
                    stderr = sys.stderr)

    # Finally - OPTIONAL - we may want to build a score database for the
    # resulting Shear Catalogs.

