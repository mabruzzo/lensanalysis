"""
This module updates extends the capability of LensTools to allow for more 
flexibility with the choice of scoring method for the while constructing the 
chi2database
"""

from __future__ import division,print_function,with_statement

from operator import mul
from functools import reduce

import numpy as np
import pandas as pd
from scipy import interpolate

from lenstools import Ensemble
from lenstools.utils.decorators import Parallelize
from lenstools.simulations.logs import logdriver
from lenstools.statistics.database import ScoreDatabase
from lenstools.statistics.constraints import _interpolate_wrapper
from emcee.ensemble import _function_wrapper

class AugmentedRbf(interpolate.Rbf):
    def __init__(self, *args, **kwargs):
        self.xi = np.asarray([np.asarray(a, dtype=np.float_).flatten()
                           for a in args[:-1]])
        self.N = self.xi.shape[-1]
        self.di = np.asarray(args[-1]).flatten()

        if not all([x.size == self.di.size for x in self.xi]):
            raise ValueError("All arrays must be equal length.")

        self.norm = kwargs.pop('norm', self._euclidean_norm)
        self.epsilon = kwargs.pop('epsilon', None)
        if self.epsilon is None:
            # default epsilon is the "the average distance between nodes" 
            # based on a bounding hypercube
            dim = self.xi.shape[0]
            ximax = np.amax(self.xi, axis=1)
            ximin = np.amin(self.xi, axis=1)
            edges = ximax-ximin
            edges = edges[np.nonzero(edges)]
            self.epsilon = np.power(np.prod(edges)/self.N, 1.0/edges.size)
        self.smooth = kwargs.pop('smooth', 0.0)

        self.function = kwargs.pop('function', 'multiquadric')
        
        b = np.zeros((len(self.di)+self.xi.shape[0]+1,))
        b[:len(self.di)] = self.di
        x = np.linalg.solve(self._A(),b)
        self.polycoef = x[len(self.di):-1][None].T
        self.offset = x[-1]
        self.nodes = x[:len(self.di)]

    def _B(self):
        r = self._call_norm(self.xi, self.xi)
        return self._init_function(r) - np.eye(self.N)*self.smooth

    def _A(self):
        B = self._B()
        tP = np.concatenate((self.xi,np.ones((self.N,))[None]),axis=0)
        zero_chunk = np.zeros((tP.shape[0],tP.shape[0]))
        A = np.block([[B,tP.T],
                      [tP,zero_chunk]])
        return A
    
    def __call__(self, *args):
        args = [np.asarray(x) for x in args]
        if not all([x.shape == y.shape for x in args for y in args]):
            raise ValueError("Array lengths must be equal")
        shp = args[0].shape
        xa = np.asarray([a.flatten() for a in args], dtype=np.float_)
        add = (xa * self.polycoef).sum(axis=0) + self.offset
        r = self._call_norm(xa, self.xi)
        return np.dot(self._function(r), self.nodes).reshape(shp) + add


def manual_rbf_interpolator_setup(emulator,function = 'multiquadric', 
                                  augmented = False, ref = None):
    used_parameters = emulator.parameter_set
    
    if '_num_bins' not in emulator._metadata:
        emulator._metadata.append("_num_bins")
    emulator._num_bins = reduce(mul, emulator.feature_set.shape[1:])

    new_shape = (emulator.feature_set.shape[0],emulator._num_bins)
    flattened_feature_set = emulator.feature_set.reshape(new_shape)
    
    # Build the interpolator
    if '_interpolator' not in emulator._metadata:
        emulator._metadata.append("_interpolator")
        
    interp = []
    for n in range(emulator._num_bins):
        args = (tuple(used_parameters.T)+(flattened_feature_set[:,n],))
        kwargs = {'function':function}
        if ref is not None:
            kwargs['epsilon']= ref._interpolator[n]().epsilon
        if augmented:
            interp.append(_interpolate_wrapper(AugmentedRbf,
                                               args=args,
                                               kwargs=kwargs))
        else:
            interp.append(_interpolate_wrapper(interpolate.Rbf,
                                               args=args,
                                               kwargs=kwargs))
    emulator._interpolator = interp


    
class BetterChi2Scorer(object):
    """
    The built in chi2 scoring method always automatically uses the
    default scoring method.
    
    This method let's you specify which scoring method you wish to use.
    """
    
    def __init__(self,function='multiquadric',augmented = False):
        self.function = function
        self.augmented = augmented
        self.features_covariance = None

    def __call__(self,sub_emulator,parameter_array,sub_feature,**kwargs):
        if "features_covariance" in kwargs:
            features_covariance = kwargs["features_covariance"]
            self.features_covariance = features_covariance
            del(kwargs["features_covariance"])
        else:
            features_covariance = self.features_covariance
            assert features_covariance is not None

        c = sub_emulator.feature_names
        sub_emulator_columns = sub_emulator[c].columns
        #sub_emulator.train(function = self.function)
        manual_rbf_interpolator_setup(sub_emulator,
                                      function = self.function, 
                                      augmented = self.augmented, 
                                      ref = None)

        temp = features_covariance[c][sub_emulator_columns].loc[c]
        sub_feature_covariance = temp.loc[sub_emulator_columns].values
        return sub_emulator.chi2(parameters=parameter_array,
                                 observed_feature=sub_feature,
                                 features_covariance=sub_feature_covariance,
                                 **kwargs)

def chi2score(emulator, parameters, data, data_covariance, nchunks, pool,
              method):

	#Score the data on each of the parameter combinations provided
	scores = emulator.score(parameters, data,
                                features_covariance=data_covariance,
                                split_chunks=nchunks,pool=pool,
                                method = method)

	#Pop the parameter columns, compute the likelihoods out of the chi2
	for p in parameters.columns:
		scores.pop(p)

        return scores,scores.apply(lambda c:np.exp(-0.5*c),axis=0)


@Parallelize.masterworker
def chi2database(db_name, parameters, specs, table_name="scores", pool=None,
                 nchunks = None, score_method = 'chi2'):
    """
    Populate an SQL database with the scores of different parameter sets with 
    respect to the data; supports multiple features.

    Parameters
    ----------
    db_name: str
        name of the database to populate
    parameters: lenstools.Ensemble
        parameter combinations to score
    specs: dict
        dictionary that should contain the emulator, data, and covariance 
        matrix of each feature to consider; each value in this dictionary must 
        be a dictionary with keys 'emulator', 'data' and 'data covariance'.
    table_name: str
        table name to populate in the database
    pool: emcee.utils.MPIPool, optional 
        MPIPool to spread the calculations over (pass None for automatic pool 
        handling)
    nchunks: int
        number of chunks to split the parameter score calculations in (one 
        chunk per processor ideally)
    score_method: str or callable
        scoring method to use (defaults to chi2). If callable, must take in the
        current instance of Emulator, the parameter array and the observed 
        feature and return a score for each parameter combination

    Notes
    -----
    This is the exact implementation provided in lenstools.statistics.database 
    except that we have added score_method.
    """
    
    #Each processor should have the same exact workload
    if nchunks is not None:
        assert not len(parameters)%nchunks
    #Database context manager
    logdriver.info(("Populating table '{0}' of score database {1}"
                    "...").format(table_name,db_name))
    with ScoreDatabase(db_name) as db:
        #Repeat the scoring for each key in the specs dictionary
        for feature_type in specs.keys():
            #Log
            logdriver.info(("Processing feature_type: {0} ({1} feature "
                            "dimensions, {2} parameter combinations)"
                            "...").format(feature_type,
                                          len(specs[feature_type]["data"]),
                                          len(parameters)))

            # score
            temp = chi2score(emulator = specs[feature_type]["emulator"],
                             parameters = parameters,
                             data = specs[feature_type]["data"],
                             data_covariance = \
                             specs[feature_type]["data_covariance"],
                             nchunks = nchunks, pool = pool,
                             method = score_method)
            chi2,likelihood = temp
            assert (chi2.columns==[feature_type]).all()

            #Add to the database
            db_chunk = parameters.copy()
            db_chunk["feature_type"] = feature_type
            db_chunk["chi2"] = chi2
            db_chunk["likelihood"] = likelihood

            db.insert(db_chunk,table_name)
