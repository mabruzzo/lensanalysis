import sys,os
import argparse,ConfigParser
import logging
import time

import lenstools as lt
from lenstools.simulations.logs import peakMemory,peakMemoryAll
import astropy.units as u

#MPI
from mpi4py import MPI
from lenstools.utils import MPIWhirlPool

"""
this script will postprocess (smooth and add noise) to convergence maps.
"""

class PathTracker(object):
    def __init__(self, noise_fname, unsmoothed_fname, smoothed_kappa_fname,
                 smoothed_noisy_fname, smoothed_noisy_catalog_fname):

        self.noise_fname = noise_fname
        self.unsmoothed_fname = unsmoothed_fname
        self.smoothed_kappa_fname = smoothed_kappa_fname
        self.smoothed_noisy_fname = smoothed_noisy_fname
        self.smoothed_noisy_catalog_fname = smoothed_noisy_catalog_fname

    def noise_fname(self,realization):
        return self.noise_fname.format(realization)

    def unsmoothed(self,realization):
        return self.unsmoothed_fname.format(realization)

    def smoothed_kappa(self,realization):
        return self.smoothed_kappa_fname.format(realization)

    def smoothed_noisy(self,realization):
        return self.smoothed_noisy_fname.format(realization)

    def smoothed_noisy_catalog(self,realization):
        return self.smoothed_noisy_catalog_fname.format(realization)

def ensure_dir_exists(paths):
    """
    Ensures that the directories exist, and if they don't, create them.
    """
    for path in paths:
        if os.path.isdir(path):
            continue
        os.makedirs(path)

def smoothing(map_id, path_tracker, overwrite, smoothing_scale):
    """
    Try to load in the smoothing scale
    """
    smoothed_kappa_fname = path_tracker.smoothed_kappa(map_id)
    if (not overwrite) and os.path.isfile(smoothed_kappa_fname):
        return lt.ConvergenceMap.load(smoothed_kappa_fname)
    else:
        conv_map = lt.ConvergenceMap.load(conv_fname)
        smoothed_kappa = conv_map.smooth(smoothing_scale)
        smoothed_kappa.save(smoothed_kappa_fname)
        return smoothed_kappa

def add_noise(map_id, path_tracker, overwrite, smoothed_kappa):
    smoothed_noisy_fname = path_tracker.smoothed_kappa(map_id)
    if (not overwrite) and os.path.isfile(smoothed_noisy_fname):
        # just return the pre-existing map
        return lt.ConvergenceMap.load(smoothed_noisy_fname)
    else:
        noise_fname = path_tracker.noise_file(map_id)
        noise = np.load(noise_fname)
        smoothed_noisy = smooth + noise
        return smoothed_noisy

def catalog_peaks(smoothed_kappa, smoothed_noisy, map_id, path_tracker,
                  overwrite=True, options={}):
    """
    Catalogs the peaks.
    One of the possible options 'both' is to catalog both peaks in 
    smoothed_kappa and smoothed_noisy. The default is False - only catalog 
    smoothed_noisy
    overwrite means nothing
    """
    catalog_both = options.get('both',default=False)
    if catalog_both:
        extent = [smoothed_kappa.data.min(), smoothed_kappa.data.max()]
        height,positions = smooth_kappa.locatePeaks(extent)
        smoothed_kappa_peaks = np.column_stack((height,positions[:,0].value,
                                                positions[:,1].value))
        
        raise NotImplementedError()
        #np.save(k_peak_cat_fname,smoothed_kappa_peaks)
    extent = [smoothed_noisy.data.min(), smoothed_noisy.data.max()]
    height,positions = smoothed_noisy.locatePeaks(extent)
    smoothed_noisy_peaks = np.column_stack((height,positions[:,0].value,
                                            positions[:,1].value))
    np.save(smoothed_noisy_catalog, smoothed_noisy_peaks)

def post_process(map_id,path_tracker,smoothing_scales,logger,overwrite=True,
                 extra_callback = [], options = {}):
    logger.info("Processing realization {:04d}".format(map_id))
    # compute smoothed_kappa
    smoothed_kappa = smoothing(map_id, path_tracker, overwrite, smoothing_scale)
    # compute smoothed_noisy
    smoothed_noisy = add_noise(map_id, path_tracker, overwrite, smoothed_kappa)

    # lets call the callback functions
    for callback in extra_callback:
        callback(smoothed_kappa, smoothed_noisy, map_id, path_tracker,
                 overwrite, options)

def dispatcher(pool, num_realizations, path_tracker, smoothing_scales, logger,
               overwrite=True, extra_callback = [], options = {}):

    if pool is None:
        first_map_realization = 0
        last_map_realization = num_realizations
        realizations_per_test = num_realizations
        fmat = (first_map_realization, last_map_realization)
        logger.debug("Post Processing realizations from {0} to {1}".format(fmat))
    else:
        message = ("Perfect load-balancing enforced, map_realizations must be a "
                   "multiple of the number of MPI tasks!")
        assert map_realizations%(pool.size+1)==0, message
        realizations_per_task = map_realizations//(pool.size+1)
	first_map_realization = realizations_per_task*pool.rank
        last_map_realization = realizations_per_task*(pool.rank+1)
        logger.debug("Task {0} will generate lensing map realizations from {1}"
                        "to {2}".format(pool.rank,first_map_realization+1,
                                        last_map_realization))
    begin = time.time()
    #Log initial memory load
    peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
    if (pool is None) or (pool.is_master()):
        logger.info(("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all "
                        "{1[1]} tasks)").format(peak_memory_task,peak_memory_all))

    for realization in range(first_map_realization,last_map_realization):
        last_timestamp = now
        post_process(map_id+1,path_tracker,smoothing_scales,logger,overwrite,
                     extra_callback = extra_callback, options = options)
        now = time.time()
		
	#Log peak memory usage to stdout
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	logger.info(("Weak lensing calculations for realization {0} completed in "
                     "{1:.3f}s").format(realization+1,now-last_timestamp))
	logger.info(("Peak memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} "
                     "tasks)").format(peak_memory_task,peak_memory_all))

    #Safety sync barrier
    if pool is not None:
        pool.comm.Barrier()

    if (pool is None) or (pool.is_master()):	
	now = time.time()
	logdriver.info("Total runtime {0:.3f}s".format(now-begin))

if __name__ == '__main__':

    # This needs to be cleaned up. Namely, we need to put some of the
    # operations inside of functions
    # we also need to do error checking with the configuration file

    
    #Parse command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",dest="verbose",action="store_true",
                        default=False,help="turn output verbosity")
    parser.add_argument("-r","--resume",dest="resume",action="store_true",
                        default = False,
                        help = ("resume previous processing - don't overwrite "
                                "existing post processed files"))
    parser.add_argument("-p","--peaks",dest="peaks",action="store_true",
                        default = False, help="turn on peak catalog creation")
    parser.add_argument("-z","--redshift", dest="z", action="store", type=float,
                        help="redshift of the source plane, required",
                        required = True)
    parser.add_argument("-c","--config",dest="config_file",action="store",
                        type=str,help="configuration file",required=True)
    parser.add_argument("id",nargs="*")

    #Parse command arguments
    cmd_args = parser.parse_args()

    # create logger
    logger = logging.getLogger("Logger")

    conh = logging.StreamHandler()
    #Verbosity level
    if cmd_args.verbose:
	logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter_console = logging.Formatter('%(levelname)s - %(message)s')
    conh.setFormatter(formatter_console)
    logger.addHandler(conh)
    
    #Initialize MPIWhirlPool
    comm = MPI.COMM_WORLD

    try:
	pool = MPIWhirlPool(comm=comm)
    except:
	pool = None
    logger.debug("Couldn't initialize MPI Pool, running in series")

    # load the configuration file

    config = ConfigParser.SafeConfigParser()
    config.read(cmd_args.config_file)

    # find allowed redshifts
    zs = map(float,zip(*config.items('Source Plane Redshifts'))[1])
    

    if cmd_args.z not in zs:
        raise ValueError("Invalid z chosen")

    # check to see if the selected redshift deviates from normal formatting
    if 'Format Deviation' in config.sections():
        problem_config_indices = map(int,
                                     zip(*config.items('Format Deviation'))[1])
        problem_zs = [zs[i] for i in problem_config_indices]
        if cmd_args.z in problem_zs:
            raise NotImplementedError("Can't handle redshifts that deviate "
                                      "from format")

    z = cmd_args.z

    # start setting up the file paths

    noise_fname= config.get('Noise Files','file_template')
    temp = '/'.join([config.get('Unsmoothed Maps', 'root_directory'),
                     config.get('Unsmoothed Maps',
                                'directory_template').format(z)])
    conv_map_fname = ''.join([config.get('Unsmoothed Maps',
                                         'file_prefix').format(z),
                              config.get('Unsmoothed Maps',
                                         'realization_portion'),
                              config.get('Unsmoothed Maps',
                                         'file_suffix')])
    unsmoothed_fname = '/'.join([temp, conv_map_fname])

    analysis_dir = config.get('Analysis Data','root_directory')
    processed_maps = '/'.join([analysis_dir,
                               'Maps',
                               config.get('Analysis Data',
                                          'directory_template').format(z)])
    smoothed_kappa_dir = '/'.join([processed_maps,'Smoothed_kappa'])
    ensure_dir_exists(smoothed_kappa_dir)
    smoothed_kappa_fname = '/'.join([smoothed_kappa_dir,conv_map_fname])
    smoothed_noisy_dir = '/'.join([processed_maps,'Smoothed_noisy'])
    ensure_dir_exists(smoothed_noisy_dir)
    smoothed_noisy_fname = '/'.join([smoothed_noisy_dir,conv_map_fname])
    source_config_dir = '/'.join([analysis_dir,'fiducial',
                                  config.get('Analysis Data',
                                             'directory_template').format(z)])
    catalog_dir = '/'.join([source_config_dir,'PeakCatalogs','KNPeaks'])
    ensure_dir_exists(catalog_dir)
    smoothed_noisy_catalog_fname = '/'.join([catalog_dir,'peaks_{04:d}r.npy'])
    

    # create an instance of path tracker
    path_tracker = PathTracker(noise_fname, unsmoothed_fname,
                               smoothed_kappa_fname, smoothed_noisy_fname,
                               smoothed_noisy_catalog_fname)

    # execute

    num_realizations = config.getint('Analysis Data','num_realizations')
    overwrite = not cmd_args.resume
    extra_callback = []
    if cmd_args.peaks:
        extra_callback.append(catalog_peaks)

    smoothing_scale = (1.0*u.arcminute)/np.sqrt(2)
