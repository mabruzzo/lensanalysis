import argparse
import ConfigParser
import os, os.path
from string import Formatter

import numpy as np
from glob2 import glob

from lenstools.statistics.ensemble import Ensemble

from lensanalysis.utility.emulator_building import  construct_emulator, \
    build_covariance_matrix, load_fiducial, build_ensembles
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
                            "Default is the minimum allowed realization to "
                            "be processed. This is one-indexed."))
parser.add_argument("--max",dest = "max_realization", action = "store",
                    default = None, type =int,
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed. This is "
                            "one-indexed."))
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
parser.add_argument("-p", dest = "p", action = "store_true", default = False,
                    help = ("Indicates that we are using MPI"))

def load_config_details(cmd_args):
    path = cmd_args.config
    assert os.path.isfile(path)

    config = ConfigParser.SafeConfigParser()
    config.read(path)
    out = {}
    num_feat_bins = config.getint("Details","num_feat_bins")
    assert num_feat_bins>0
    out['num_feat_bins'] = num_feat_bins

    if config.has_option("Details","combine_neighboring_bins"):
        combine_neighboring_bins = config.getint("Details",
                                                 "combine_neighboring_bins")
        if combine_neighboring_bins <0:
            raise ValueError()
        elif combine_neighboring_bins>1:
            assert num_feat_bins % combine_neighboring_bins ==0
        out['combine_neighboring_bins'] = combine_neighboring_bins
    else:
        out['combine_neighboring_bins'] = 0
    
    num_tomo_bins = config.getint("Details","num_tomo_bins")
    assert num_feat_bins>0
    out['num_tomo_bins'] = num_tomo_bins

    out['cross_statistic'] = config.getboolean("Details","cross_statistic")

    if out['cross_statistic']:
        temp = config.get("Details", "specific_cross_bins")
        l = [int(elem) for elem in temp.split(',')]
        if len(l) == 0:
            raise ValueError("Must specify at least one value for "
                             "specific_cross_bins")
        elif len(l) == 1 and l[0] == -1:
            out["specific_cross_bins"] = None
        else:
            for elem in l:
                assert elem>=0
            out["specific_cross_bins"] = l
    else:
        out["specific_cross_bins"] = None

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

    for par_min,par_max in [("Om_min","Om_max"), ("w_min","w_max"),
                            ("sigma8_min","sigma8_max")]:
        out[par_min] = config.getfloat("Interpolation",par_min)
        out[par_max] = config.getfloat("Interpolation",par_max)
        if out[par_min] >= out[par_max]:
            if out[par_min] == out[par_max] and par_min == "w_min":
                continue
            raise ValueError("{:s} must be less than {:s}".format(par_min,
                                                                  par_max))

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
        assert pca <= config_dict['num_feat_bins']*config_dict['num_tomo_bins']
        if len(tomo_bins) > 1:
            assert pca <= config_dict['num_feat_bins']
    else:
        assert pca <= config_dict['num_feat_bins']
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


    num_feat_bins = config_dict["num_feat_bins"]
    combine_neighbor = config_dict['combine_neighboring_bins']
    if combine_neighbor >1:
        # sanity check
        assert num_feat_bins %combine_neighbor == 0
        num_feat_bins = num_feat_bins //combine_neighbor

    for tomo_bin in tomo_bins:
        if tomo_bin == -1:
            feature_name = "TomoPeaks"
            feature_index = []
            bin_num = None

            if config_dict["cross_statistic"]:
                if config_dict["specific_cross_bins"] is not None:
                    bin_num = np.array(config_dict["specific_cross_bins"]) + 1
                    # sorting is important!
                    bin_num.sort()
                    count = 0

                for i in range(1,config_dict["num_tomo_bins"]+1):
                    for j in range(i,config_dict["num_tomo_bins"]+1):

                        if bin_num is not None:
                            count += 1
                            if count not in bin_num:
                                continue

                        if i == j:
                            for k in range(num_feat_bins):
                                feature_index.append("z{:d}_{:d}".format(i,k))
                        else:
                            for k in range(num_feat_bins):
                                feature_index.append(("z{:d},{:d}_{:d}"
                                                      ).format(i,j,k))
            else:
                for i in range(1,config_dict["num_tomo_bins"]+1):
                    for j in range(num_feat_bins):
                        feature_index.append("z{:d}_{:d}".format(i,j))
            

        else:
            feature_name = "z{:02d}".format(tomo_bin)
            feature_index = ["z{:d}_{:d}".format(tomo_bin,j) \
                             for j in range(num_feat_bins)]
            if config_dict["cross_statistic"]:
                # first convert tomo_bin to zero_indexed
                ind = tomo_bin - 1
                bin_num = (ind*(2*config_dict["num_tomo_bins"] + 1 - ind)//2
                           + ind - ind)
                # finally conver bin_num from zero-indexing back to one-indexing
                bin_num-=1
            else:
                bin_num = tomo_bin
        # build the emulator
        temp = construct_emulator(config_dict['paths'], parameter_index,
                                  feature_index = feature_index,
                                  fname_prefix=fname_prefix,
                                  bin_num=bin_num, pca = pca_bins,
                                  pca_scale = None,
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


# Below I have included some functions for identifying the parameter values
# during MPI processing

# this needs to be tested slightly better and the leftmost value should
# probably be improved

def prepare_col(ind_start,vals,length,cols_to_right=0,leftmost=False):
    if cols_to_right:
        vals = np.repeat(vals,cols_to_right*len(vals))
    l_vals = len(vals)
    if l_vals>ind_start:
        roll = l_vals - ind_start
    elif l_vals<ind_start:
        roll = l_vals - (ind_start % l_vals)
    else:
        roll = 0
    cur = np.roll(vals,roll)
    
    if leftmost:
        return cur[:length]
    num_tiles = length // l_vals
    if length % l_vals:
        num_tiles +=1
    return np.tile(cur,num_tiles)[:length]

def prepare_last_col(ind_start,si_8_vals,length):
    return prepare_col(ind_start,si_8_vals,length)

def prepare_middle_col(ind_start,w_vals,length):
    return prepare_col(ind_start,w_vals,length,1)

def prepare_first_col(ind_start,om_vals,length):
    return prepare_col(ind_start,om_vals,length,2,True)



def load_realization_ind(cmd_args, config_dict):
    """
    Load in the details about individual realizations.
    """
    if not cmd_args.indiv_realizations:
        return [None]

    min_r, max_r = cmd_args.min_realization, cmd_args.max_realization
    if min_r is None:
        min_r = 1
    if max_r is None:
        # let's load in the fiducial temporarily - we just want to see how many
        # realizations there are. We will properly load it in again later. 
        fid = build_ensembles([config_dict['fiducial_obs_path']])[1][0]
        max_r = fid.values.shape[0]

    assert isinstance(min_r,int) and isinstance(max_r,int)
    assert 1 <= min_r <=max_r
    return range(min_r,max_r+1)


def driver(cmd_args):
    config_dict = load_config_details(cmd_args)
    tomo_bins = identify_tomo_bins(cmd_args,config_dict)
    
    pca_bins = determine_pca_bins(cmd_args,config_dict,tomo_bins)
    # for now, we will not worry about the following 2 parameters
    realization_ind = load_realization_ind(cmd_args,config_dict)
    db_name_template = get_save_file(cmd_args,realization_ind)
    
    for elem in realization_ind:
        specs,ensemble_params,param_name = prepare_specs(tomo_bins,config_dict,
                                                         pca_bins = pca_bins,
                                                         realization_ind = elem)

        if cmd_args.p:
            assert len(realization_ind) == 1
            assert len(specs) == 1
            raise NotImplementedError()
        else:
            
            # probably should make the following adjustable
            # create the grid on which we will interpolate
            num_axis = config_dict['num_samples']
            if config_dict["w_min"] != config_dict["w_max"]:
                p = np.array(np.meshgrid(np.linspace(config_dict["Om_min"],
                                                     config_dict["Om_max"],
                                                     num_axis),
                                         np.linspace(config_dict["w_min"],
                                                     config_dict["w_max"],
                                                     num_axis),
                                         np.linspace(config_dict["sigma8_min"],
                                                     config_dict["sigma8_max"],
                                                     num_axis),
                                         indexing="ij")).reshape(3,
                                                                 num_axis**3).T
                test_parameters = Ensemble(p,columns=param_name)
            else:
                param_dict = {"Om" : np.linspace(config_dict["Om_min"],
                                                 config_dict["Om_max"],
                                                 num_axis),
                              "sigma8" : np.linspace(config_dict["sigma8_min"],
                                                     config_dict["sigma8_max"],
                                                     num_axis)}
                test_parameters = Ensemble.meshgrid(param_dict,
                                                    {"Om" : 0,
                                                     "sigma8": 1})
                test_parameters["w"] = -1.0
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
                raise NotImplementedError("Have not defined behavior for when "
                                          "the database already exists")
            chi2database(db_name, test_parameters, specs, table_name="scores",
                         pool=pool, nchunks=nchunks,
                         score_method = config_dict['method'])

if __name__ == '__main__':
    cmd_args = parser.parse_args()
    driver(cmd_args)
