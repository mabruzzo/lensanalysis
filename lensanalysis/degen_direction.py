from math import sqrt
import numpy as np
from scipy.optimize import minimize_scalar

"""
Here we implement functions to find the alpha that minimizes the scatter in 
Sigma_8 = sigma_8 * (Omega_m / C)^alpha.
"""

def _compute_score(alpha, sigma_8, red_omega, likelihood):
    """
    Computes the ratio of sqrt(V)/E for Sigma_8 which is given by 
    Sigma_8 = sigma_8 * (Omega_m/const)^alpha

    Here red_omega is just an Omega_m/constant

    """
    val = sigma_8 * np.power(red_omega,alpha)
    expectation = np.average(val,weights=likelihood)
    v = np.average(np.square(val-expectation),weights=likelihood)
    return sqrt(v)/expectation

def compute_alpha_Petri(sigma_8,omega_m,likelihood, const, method = 'brent',
                        tol=0.001, maxiter = 500, bounds = None,
                        full_result = False, options=None):
    """
    Determines the value of alpha for Sigma_8 = sigma_8 * (Omega_m/const)^alpha

    The method employed is described in Petri et al. 2015

    Parameters
    ----------
    sigma_8 : array_like
        The values of sigma_8 that the likelihood grid is sampled at.
    omega_m : array_like
        The values of Omega_m that the likelihood grid is sampled at. Must be 
        the same shape as sigma_8.
    likelihood : array_like
        The likelihood values of various likelihood combinations sampled on a 
        grid.
    const : float
        The constant used in the equation of Sigma_8
    method : string or callable
        The name of the scipy optimization method being used or a custom 
        minimization function. Options are 'Brent', 'Bounded', and 'Golden'. 
        Alternatively, supply a callable object.
    tol : float,optional
        The tolerance for convergence. For 'brent' and 'golden' this is the 
        relative error in xopt for convergence while for 'bounded' it is the 
        absolute error in xopt needed for convergence.
    maxiter : int,optional
        Maximum number of iterations
    bounds : sequence,optional
        The minimum and maximum value to use to constrain alpha with the 
        'bounded' method.
    full_result : bool,optional
        Whether or not to return the direct output of the minimization function.
        This is always done when a callable method is passed in. If this is 
        False, and the minimization is unsuccessful, None is returned.
    options: dict, optional
        Optional keyword arguments to pass a callable method.
    """

    if callable(method):
        meth = 'custom'
    else:
        meth = method

    # need to check shapes of sigma_8, omega_m, likelihood
    if len(likelihood.shape)>2:
        raise NotImplementedError("Not currently prepared for likelihoods "
                                  "sampled in more than 2 dimensions")

    assert sigma_8.shape == omega_m.shape
    assert likelihood.shape == sigma_8.shape

    if len(sigma_8.shape)>1:
        sigma_8 = sigma_8.flatten()
        omega_m = omega_m.flatten()
        likelihood = likelihood.flatten()
    # deal with NaNs
    w = np.logical_not(np.isnan(likelihood))
    sigma_8 = sigma_8[w]
    omega_m = omega_m[w]
    likelihood = likelihood[w]

    # need to compute reduced omega_m - just omega divided by constant
    # this will save us time by not needing to recalculate omega_m/const
    red_omega = omega_m/float(const)

    # now use a scipy optimization method or a user suplied method to perform
    # optimization


    if meth == 'custom':
        res = minimize_scalar(_compute_score, method = method, bracket = None,
                              bounds = bounds,
                              args = (sigma_8, red_omega, likelihood),
                              tol = tol,options = options)
        return res
    else:
        options = {'maxiter' : maxiter}
        res = minimize_scalar(_compute_score, method = method, bracket = None,
                              bounds = bounds,
                              args = (sigma_8, red_omega, likelihood),
                              tol = tol,options = options)
    if full_result:
        return res
    elif res.success:
        return result.x
    else:
        return None
