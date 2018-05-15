import argparse
import ConfigParser
import os, os.path
from string import Formatter

import numpy as np
from glob2 import glob

from lenstools.statistics.ensemble import Ensemble

from lensanalysis.utility.emulator_building import  construct_emulator, \
    build_covariance_matrix, load_fiducial
from lensanalysis.utility.scoring import BetterChi2Scorer, chi2database



"""
Use the chi2database to construct a database of the chi-squared values for 
different parameter combinations.
"""

parser = argparse.ArgumentParser(description = "Construct a score database")


parser.add_argument('-o', '--observed', dest = 'observed', action = 'store',
                    default = None,
                    help = ('the path to the file with all of the the average '
                            'peak counts for the observed features. If this is '
                            'not provided, then we use the fiducial_obs_path'
                            'supplied in the configuration file.'))
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
parser.add_argument("--config", dest = "config",action = "store",
                    help = ("Configuration file that include details about "
                            "what data should be used to score the "
                            "observations."))
parser.add_argument("--pca", dest = "pca", action = "store",
                    default = None, type = int,  
                    help=("The number of bins to use for PCA. If not "
                          "specified, then PCA is not performed"))
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
                            "values for --indiv_tomo_bin must be specified."))
parser.add_argument("--indiv_tomo_bin", dest = "indiv_tomo_bin",
                    action = "store", type = int, default = None, nargs = '*',
                    help = ("Specify this option only if you want to build a "
                            "score databases for individual tomographic bins. "
                            "This expects integers refering to the index (one-"
                            "indexed) of tomographic bins that have been "
                            "ordered by increasing redshift. If not specified, "
                            "the score database is only built for the "
                            "full tomographic dataset. If the option is "
                            "specified, but no bins are given, then all "
                            "indivdual tomographic bins are processed."))

def load_config_details(cmd_args):
    path = cmd_args.config
    assert os.path.isfile(path)

    config = ConfigParser.SafeConfigParser()
    config.read(path)
    out = {}
    num_peak_bins = config.getint("Details","num_peak_bins")
    assert num_peak_bins>0
    out['num_peak_bins'] = num_peak_bins

    if config.has_option("Details","combine_neighboring_bins"):
        combine_neighboring_bins = config.getint("Details",
                                                 "combine_neighboring_bins")
        if combine_neighboring_bins <0:
            raise ValueError()
        elif combine_neighboring_bins>1:
            assert num_peak_bins % combine_neighboring_bins ==0
        out['combine_neighboring_bins'] = combine_neighboring_bins
    else:
        out['combine_neighboring_bins'] = 0
    
    num_tomo_bins = config.getint("Details","num_tomo_bins")
    assert num_peak_bins>0
    out['num_tomo_bins'] = num_tomo_bins

    # load in the covariance matrix array
    covariance_cosmo_path = config.get("Details","covariance_cosmo_path")
    assert os.path.isfile(covariance_cosmo_path)
    out['covariance_cosmo_path'] = covariance_cosmo_path

    # load in the observed feature
    if cmd_args.observed is None:
        str_insert = "The fiducial_obs_path configuration option "
        observed = config.get("Details","fiducial_obs_path")
    else:
        str_insert = "The -o/--observed option "
        observed = cmd_args.observed

    if not os.path.isfile(observed):
        raise ValueError("{:s}must point to an existing file.".format(observed))
    out['fiducial_obs_path'] = observed

    # get the directory where the average peak counts of the various sampled
    # cosmologies are saved.
    root_dir = config.get("Details","root")
    assert os.path.isdir(root_dir)
    paths = sorted(glob('{}/*.npy'.format(root_dir)))
    assert len(paths)>0
    out['paths'] = paths

    # finally get the number of characters preceeding the cosmo information
    # in the sampled cosmology average peak count files.
    root_fname_prefix = config.getint("Details","root_fname_prefix")
    assert root_fname_prefix > 0
    out['root_fname_prefix'] = root_fname_prefix


    # now we will just extract some details related to the interpolation
    # first lets get the number of samples per parameter
    out['num_samples'] = config.getint("Interpolation","num_samples")
    assert out['num_samples'] > 0

    # next lets configure Chi2Scoring method
    augmented = config.getboolean("Interpolation","augmented_Rbf")
    basis_function = config.get("Interpolation", "basis_function")
    assert basis_function in ['multiquadric', 'inverse', 'gaussian', 'linear',
                              'cubic', 'quintic', 'thin_plate']
    out['method'] = BetterChi2Scorer(function = basis_function,
                                     augmented = augmented)

    return out

def identify_tomo_bins(cmd_args, config_dict):
    temp = cmd_args.indiv_tomo_bin
    if temp is None:
        # we do not want to process individual tomo bins
        tomo_bins = []
    elif len(temp) == 0:
        # we want to process all individual tomo bins
        tomo_bins = range(1, config_dict['num_tomo_bins']+1)
    else:
        tomo_bins = temp

    for elem in tomo_bins:
        assert elem >0

    if not cmd_args.no_tomo:
        tomo_bins.append(-1)

    if len(tomo_bins) == 0:
        raise ValueError("At least one individual tomographic bin must be "
                         "specified or we must process all tomographic bins.")
    return tomo_bins

def determine_pca_bins(cmd_args,config_dict,tomo_bins):
    pca = cmd_args.pca
    if pca is None:
        return None
    assert pca > 0
    if -1 in tomo_bins:
        assert pca <= config_dict['num_peak_bins']*config_dict['num_tomo_bins']
        if len(tomo_bins) > 1:
            assert pca <= config_dict['num_peak_bins']
    else:
        assert pca <= config_dict['num_peak_bins']
    return pca

def get_save_file(cmd_args,realizations):
    template = cmd_args.save

    iterable = Formatter().parse(template)

    num = 0
    for _, field_name, format_spec, _ in iterable:
        if field_name is not None:
            num +=1

    if len(realizations) >1:
        if num != 1:
            raise ValueError("For individual realizations, the save_file must "
                             "have 1 formatter template.")
    elif len(realizations) == 1:
        if realizations[0] == None:
            if num != 0:
                raise ValueError("For a model without any realizations, there "
                                 "must not be any realizations")
        elif num !=1:
            raise ValueError("For individual realizations, the save_file must "
                             "have 1 formatter template.")
    return template

def prepare_specs(tomo_bins,config_dict, pca_bins = None,
                  realization_ind = None):
    """
    This prepares the specifications to be processed.
    """
    specs = {}
    # note a tomo_bin value of -1 indicates that we want full tomography
    parameter_index = ["Om","w","sigma8"]
    fname_prefix=config_dict['root_fname_prefix']
    cov_path = config_dict['covariance_cosmo_path']
    fid_path = config_dict['fiducial_obs_path']


    num_peak_bins = config_dict["num_peak_bins"]
    combine_neighbor = config_dict['combine_neighboring_bins']
    if combine_neighbor >1:
        # sanity check
        assert num_peak_bins %combine_neighbor == 0
        num_peak_bins = num_peak_bins //combine_neighbor

    for tomo_bin in tomo_bins:
        if tomo_bin == -1:
            feature_name = "TomoPeaks"
            feature_index = []
            for i in range(1,config_dict["num_tomo_bins"]+1):
                for j in range(num_peak_bins):
                    feature_index.append("z{:d}_{:d}".format(i,j))
            bin_num = None
            
        else:
            feature_name = "z{:02d}".format(tomo_bin)
            feature_index = ["z{:d}_{:d}".format(tomo_bin,j) \
                             for j in range(num_peak_bins)]
            bin_num = tomo_bin
        # build the emulator
        temp = construct_emulator(config_dict['paths'], parameter_index,
                                  feature_index = feature_index,
                                  fname_prefix=fname_prefix,
                                  bin_num=bin_num, pca = pca_bins,
                                  feature_name = feature_name,
                                  combine_neighbor = combine_neighbor)
        emulator, ensemble_params,pca_basis = temp

        # load in the covariance matrix
        emulator_col = emulator[[feature_name]].columns
        features_covariance = build_covariance_matrix(cov_path,
                                                      emulator_col,
                                                      fname_prefix=0,
                                                      bin_num=bin_num,
                                                      pca_basis = pca_basis,
                                                      pca = pca_bins,
                                                      combine_neighbor = \
                                                      combine_neighbor)
        # load in the observation
        fiducial_obs = load_fiducial(fid_path, emulator_col, fname_prefix=0,
                                     bin_num = bin_num, pca_basis = pca_basis,
                                     pca = pca_bins,
                                     indiv_realization = realization_ind,
                                     combine_neighbor = combine_neighbor)
            
        specs[feature_name] = {"emulator" : emulator,
                               "data" : fiducial_obs,
                               "data_covariance" : features_covariance}
    return specs, ensemble_params, emulator.parameter_names

def driver(cmd_args):
    config_dict = load_config_details(cmd_args)
    tomo_bins = identify_tomo_bins(cmd_args,config_dict)
    
    pca_bins = determine_pca_bins(cmd_args,config_dict,tomo_bins)
    # for now, we will not worry about the following 2 parameters
    realization_ind = [None]
    db_name_template = get_save_file(cmd_args,realization_ind)

    for elem in realization_ind:
        specs,ensemble_params,param_name = prepare_specs(tomo_bins,config_dict,
                                                         pca_bins = pca_bins,
                                                         realization_ind = elem)

        # probably should make the following adjustable
        # create the grid on which we will interpolate - probably should make
        # this customizable
        num_axis = config_dict['num_samples']

        p = np.array(np.meshgrid(np.linspace(0.2,0.5,num_axis),
                                 np.linspace(-1.5,-0.5,num_axis),
                                 np.linspace(0.5,1.2,num_axis),
                                 indexing="ij")).reshape(3,num_axis**3).T
        test_parameters = Ensemble(p,columns=param_name)
        nchunks = 1
        pool = None
        nchunuks = None

        if elem is None:
            db_name = db_name_template
        else:
            db_name = db_name_template.format(elem)
        print db_name
        if os.path.isfile(db_name):
            # if we do not handle this, then the new scores will simply be
            # added to the old scores
            raise NotImplementedError("Have not defined behavior for when the "
                                      "emulator already exists")
        chi2database(db_name, test_parameters, specs, table_name="scores",
                     pool=pool, nchunks=nchunks,
                     score_method = config_dict['method'])

if __name__ == '__main__':
    cmd_args = parser.parse_args()
    driver(cmd_args)
