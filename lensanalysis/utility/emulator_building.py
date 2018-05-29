import numpy as np
import pandas as pd

from lenstools.statistics.constraints import Emulator
from lenstools.statistics.ensemble import Ensemble,Series
from lenstools.utils.algorithms import pca_transform

def build_ensembles(paths,prefix=0,bin_num = None,
                    combine_neighbor = 0):
    """
    bin_num is one-indexed!
    """
    ensembles = []
    cosmo_names = []
    for i,path in enumerate(paths):
        cosmo_name = path.split('/')[-1][prefix:-4]
        temp = np.load(path)
        if combine_neighbor>=1:
            assert(temp.shape[-1] % combine_neighbor == 0)
            
            new_shape = [elem for elem in temp.shape[:-1]]
            new_shape.append(temp.shape[-1]//combine_neighbor)
            new_shape.append(combine_neighbor)

            temp = np.sum(temp.reshape(new_shape),axis=-1)

        
        if len(temp.shape)>1:
            if bin_num is not None and not isinstance(bin_num,np.ndarray):
                e = Ensemble(temp[bin_num-1,...])

            else:
                
                if bin_num is not None:                    
                    cur_temp = temp[(bin_num-1,)]
                else:
                    cur_temp = temp

                temp_l = []
                for i in range(cur_temp.shape[0]):
                    temp_l.append(cur_temp[i,...])
                # Simply reshaping does not work for 3D input
                e = Ensemble(np.hstack(temp_l))
                #print temp.shape,cur_temp.shape,e.values.shape
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
                       fname_prefix=None, bin_num=None, pca = None,
                       pca_scale = None, feature_name = 'features',
                       combine_neighbor = 0):
    """
    Loads the features from the sampled cosmologies from a numpy array and then 
    constructs the emulator.
    """
    cosmo_names, ensembles = build_ensembles(feature_paths,
                                             prefix = fname_prefix,
                                             bin_num = bin_num,
                                             combine_neighbor =combine_neighbor)
    om_m, om_l, w, si = parse_parameter_vals(cosmo_names)
    temp = np.column_stack((om_m,om_l))
    assert (np.sum(temp,axis=1)==1).all()
    sampling_parameters = np.column_stack((om_m,w,si))
    ensemble_params = Ensemble(sampling_parameters,
                               columns = parameter_index)
    features = np.zeros((len(ensembles),ensembles[0].shape[0]))
    for i,ensemble in enumerate(ensembles):
        features[i,:] = ensemble.as_matrix()[:,0]

    if pca is not None:
        # it turns out that the mutliplication by
        # 1/np.sqrt(features.shape[0]-1.0) is part of the steps towards PCA
        feat = Ensemble(features)
        if pca_scale is None:
            scale = None
        elif pca_scale == 'mean':
            scale = feat.mean(0)
        elif pca_scale == 'unscaled':
            scale = 1.0
        else:
            raise ValueError("Invalid pca_scale keyword")

        # I cannot currently get the pca_basis to work properly
        pca_basis = feat.principalComponents(scale=scale)
        features = pca_transform(features,pca_basis,pca).values
        feature_index = None
    else:
        pca_basis = None

    # now we construct the feature index
    if feature_index is None:
        feature_index = range(features.shape[1])

    feature_index = Series.make_index(pd.Index(feature_index,
                                               name=feature_name))

    #print features
    #print ensemble_params
    #print feature_index
    #print parameter_index
    #print pca_basis
    emulator = Emulator.from_features(features,ensemble_params.as_matrix(),
                                      parameter_index = parameter_index,
                                      feature_index = feature_index)
    emulator.train()
    return emulator, ensemble_params, pca_basis

def build_covariance_matrix(path, emulator_columns, fname_prefix=None,
                            bin_num=None, pca_basis = None, pca = None,
                            estimation_correction = True, combine_neighbor = 0):
    fiducial = build_ensembles([path],prefix=fname_prefix,bin_num=bin_num,
                               combine_neighbor = combine_neighbor)[1][0]
    if pca_basis is not None or pca is not None:
        if pca_basis is not None and pca is not None:
            fiducial = Ensemble(pca_transform(fiducial.values.astype('float64'),
                                              pca_basis,pca).values)
        else:
            raise ValueError()
    

    """
    Now we will correct the covariance matrix for the fact that it was 
    estimated using a finite number of realizations.

    This correction is important for when we invert the covariance matrix.
    We are doing it here instead of after our inversion because we are unable 
    to do so after the inversion.

    This only works right now because we only ever use the covariance matrix 
    during this inversion operation. If that changes, I don't know if we should 
    still be doing this.


    Anyway, in many papers (e.g. Zorrilla Matilla et al. 2016) and in LensTools
    the inverse covariance matrix, Cinv is rescaled:
    out = (N-d-2)/(N-1) * Cinv
    where out is the corrected inverse covariance matrix, N is the number of 
    realizations, and d is the dimensions of the summary statistic. (I have 
    double checked that this is the operation performed by LensTools).

    I have confirmed that multiply a scalar by Cinv is equivalent to dividing 
    the covariance matrix by the same scalar before inversion.
    Thus we will do:
    Cov_out = (N-1)/(N-d-2)*Cov
    """
    if estimation_correction:
        nr, d = float(fiducial.values.shape[0]), float(fiducial.values.shape[1])
        # I have confirmed this updates values in place
        temp = (nr-1.)/(nr-d-2.)
        features_covariance = Ensemble(fiducial.covariance().as_matrix()*temp,
                                       index=emulator_columns,
                                       columns=emulator_columns)
    else:
        features_covariance = Ensemble(fiducial.covariance().as_matrix(),
                                       index=emulator_columns,
                                       columns=emulator_columns)
    return features_covariance

def load_fiducial(path, emulator_columns, fname_prefix=None, bin_num=None,
                  pca_basis = None, pca = None, indiv_realization = None,
                  combine_neighbor = 0):
    """
    Loads the fiducial cosmology from a numpy array and constructs the 
    "observed features" using the correct column names to match the emulator.
    """
    fiducial = build_ensembles([path],prefix=fname_prefix,bin_num=bin_num,
                               combine_neighbor = combine_neighbor)[1][0]
    if indiv_realization is not None:
        raise NotImplementedError("Not currently capable of handling "
                                  "individual realizations")
    if pca_basis is not None or pca is not None:
        if pca_basis is not None and pca is not None:
            fiducial = Ensemble(pca_transform(fiducial.values.astype('float64'),
                                              pca_basis,pca).values)
        else:
            raise ValueError()
    fiducial_obs = Series(fiducial.mean(axis=0).as_matrix(),
                          index=emulator_columns)
    return fiducial_obs
