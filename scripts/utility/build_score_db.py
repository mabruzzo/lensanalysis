import argparse

import numpy as np
import matplotlib.pyplot as plt
import lenstools as lt
from lenstools.statistics.constraints import Emulator
from lenstools.statistics.ensemble import Ensemble,Series
from lenstools.statistics.database import chi2database
from glob2 import glob
import pandas as pd

from mpi4py import MPI
from lenstools.utils import MPIWhirlPool

def build_ensembles(paths,prefix=0,bin_num = None,combine_neighbor = False):
    """
    bin_num is one-indexed!
    """
    ensembles = []
    cosmo_names = []
    for i,path in enumerate(paths):
        cosmo_name = path.split('/')[-1][prefix:-4]
        temp = np.load(path)
        if combine_neighbor:
            assert(temp.shape[-1] % 2 == 0)
            temp = temp[...,::2]+temp[...,1::2]
        if bin_num is not None and len(temp.shape)>1:
            e = Ensemble(temp[bin_num-1,:])
        elif len(temp.shape)>1:
            temp_l = []
            for i in range(temp.shape[0]):
                temp_l.append(temp[i,...])
            # Simply reshaping does not work for 3D input
            e = Ensemble(np.hstack(temp_l))
        else:
            assert bin_num is None or bin_num == 1
            e = Ensemble(temp)
        ensembles.append(e)
        cosmo_names.append(cosmo_name)
    return cosmo_names,ensembles


def parse_parameter_vals(names):
    omega_ms = []
    omega_ls = []
    sigma_8s = []
    ws = []
    for name in names:
        first, second, third, last = name.split('_')
        omega_ms.append(float(first[2:]))
        omega_ls.append(float(second[2:]))
        ws.append(float(third[1:]))
        sigma_8s.append(float(last[2:]))
    return omega_ms, omega_ls, ws, sigma_8s


def construct_emulator(feature_paths, parameter_index, feature_index = None,
                       fname_prefix=None, bin_num=None):
    """
    Loads the features from the sampled cosmologies from a numpy array and then 
    constructs the emulator.
    """
    cosmo_names, ensembles = build_ensembles(paths,prefix=fname_prefix,
                                             bin_num = bin_num)
    om_m, om_l, w, si = parse_parameter_vals(cosmo_names)
    temp = np.column_stack((om_m,om_l))
    assert (np.sum(temp,axis=1)==1).all()
    sampling_parameters = np.column_stack((om_m,w,si))
    ensemble_params = lt.Ensemble(sampling_parameters,
                                  columns = parameter_index)
    features = np.zeros((len(ensembles),ensembles[0].shape[0]))
    for i,ensemble in enumerate(ensembles):
        features[i,:] = ensemble.as_matrix()[:,0]

    emulator = Emulator.from_features(features,ensemble_params.as_matrix(),
                                      parameter_index = parameter_index,
                                      feature_index = feature_index)
    emulator.train()
    return emulator, ensemble_params

def load_fiducial(path, emulator_columns, fname_prefix=None, bin_num=None):
    """
    Loads the fiducial cosmology from a numpy array, constructs the "observed 
    features" and the covariance matrix using the correct column names to match
    the emulator.
    """
    fiducial = build_ensembles([path],prefix=fname_prefix,bin_num=bin_num)[1][0]
    fiducial_obs = Series(fiducial.mean(axis=0).as_matrix(),
                          index=emulator_columns)
    features_covariance = Ensemble(fiducial.covariance().as_matrix(),
                                   index=fiducial_obs.index,
                                   columns=fiducial_obs.index)
    return fiducial_obs,features_covariance

"""
Use the chi2database to construct a database of the chi-squared values for 
different parameter combinations.
"""

parser = argparse.ArgumentParser(description = "Construct a score database")

parser.add_argument('-r', '--root', dest = 'root', action = 'store',
                    help = ("the root directory where files with the average" 
                            "peak counts for all sampled cosmologies are "
                            "stored."))
parser.add_argument('-o', '--observed', dest = 'observed', action = 'store',
                    help = ('the path to the file with all of the the average '
                            'peak counts for the observed features. For now, '
                            'the feature covariance matrix must be built from '
                            'this.'))
parser.add_argument("--min",dest = "min_realization", action = "store",
                    default = None, type = int,
                    help = ("Specify the minimum realization to process."
                            "Default is the maximum allowed realization to "
                            "be processed."))
parser.add_argument("--max",dest = "max_realization", action = "store",
                    default = None, type =int,
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed."))
parser.add_argument("-s", "--save", dest = "save", action="store",
                    help = ("Specify the path of the file where the resulting "
                            "database will be saved. It is expected to have a "
                            "formatting character if being applied to "
                            "individual realizations."))
parser.add_argument("--pca", dest = "pca", action = "store",
                    default = None, type = int,
                    help=("The number of bins to use for PCA. If not "
                          "specified, then PCA is not performed"))
parser.add_argument("--cov", dest = "cov", action = "store", default = None,
                    help = ("File that contains the covariance matrix."))
parser.add_argument("--indiv_realizations", dest = "indiv_realizations",
                    action = "store_true", default = False,
                    help = ("Indicates whether or not we should compute the "
                            "score database for individual realizations "
                            "independently."))
parser.add_argument("--no_tomo", dest = "no_tomo", action = "store_true",
                    default = False,
                    help = ("This option explicitly indicates that the full "
                            "tomographic dataset should not be included in the "
                            "resulting database. If this is specified then "
                            "values for --indiv_tomo_bin must be specified.")
parser.add_argument("--indiv_tomo_bin", dest = "indiv_tomo_bin",
                    action = "store", type = int, default = None,
                    help = ("Specify this option only if you want to build a "
                            "score databases for individual tomographic bins. "
                            "This expects integers refering to the index (one-"
                            "indexed) of tomographic bins that have been "
                            "ordered by increasing redshift. If not specified, "
                            "the score database is only built for the "
                            "full tomographic dataset."))


if __name__ == '__main__':
    num_peak_bins = 30
    num_tomo_bins = 5
    feature_index = []
    for i in range(1,num_tomo_bins+1):
        for j in range(num_peak_bins):
            feature_index.append("z{:d}_{:d}".format(i,j))
    parameter_index = ["Om","w","sigma8"]

    root = '/home/mabruzzo/papers/photoz_wl/data/Replicating/30_bin/sampled'
    paths = sorted(glob('{}/*.npy'.format(root)))
    fid_path = ("/home/mabruzzo/papers/photoz_wl/data/Replicating/30_bin/"
                "fiducial_cosmo/30_bin_Shear.npy")

    photoz_fid_path = ("/home/mabruzzo/papers/photoz_wl/data/Replicating/"
                       "30_bin/fiducial_cosmo/30_bin_constPosBias_ppz.npy")

    # construct and train the emulator
    emulator, ensemble_params = construct_emulator(paths, parameter_index,
                                                   feature_index =feature_index,
                                                   fname_prefix=7, bin_num=None)
    # get the observed matrix and the covariance matrix
    emulator_col = emulator[["features"]].columns
    fiducial_obs,features_covariance = load_fiducial(fid_path, emulator_col,
                                                     fname_prefix=7)

    fiducial_obs = load_fiducial(photoz_fid_path, emulator_col,
                                 fname_prefix=7)[0]

    # create the grid on which we will interpolate
    num_axis = 150

    p = np.array(np.meshgrid(np.linspace(0.2,0.5,num_axis),
                             np.linspace(-1.5,-0.5,num_axis),
                             np.linspace(0.5,1.2,num_axis),
                             indexing="ij")).reshape(3,num_axis**3).T
    test_parameters = Ensemble(p,columns=emulator.parameter_names)

    db_name = ("/home/mabruzzo/papers/photoz_wl/data/Replicating/30_bin/"
               "score_db/tomo_peak_constPosBias_ppz_150_samples.sqlite")

    specs = {"features": {"emulator" : emulator,
                          "data" : fiducial_obs,
                          "data_covariance" : features_covariance}}

    #Initialize MPIWhirlPool
    #comm = MPI.COMM_WORLD
    #pool = MPIWhirlPool(comm=comm)
    nchunks = 1
    pool = None
    nchunuks = None
    chi2database(db_name,test_parameters,specs,table_name="scores",pool=pool,
                 nchunks=nchunks)
