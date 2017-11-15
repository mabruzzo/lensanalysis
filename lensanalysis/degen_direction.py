from math import sqrt
import numpy as np
from scipy.optimize import minimize_scalar

"""
Here we implement functions to find the alpha that minimizes the scatter in 
Sigma_8 = sigma_8 * (Omega_m / C)^alpha.
"""


def weighted_mean(vals,likelihood_vals):
    """
    Computes the weighted mean of parameter_grid using likelihood_grid
    """
    return np.average(vals,weights=likelihood_vals)


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

def compute_alpha_Petri(sigma_8,omega_m,likelihood, const, method, tol,
                        maxiter = 500):
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
    method : string
        The name of the scipy optimization method being used
    tol : float
        The tolerance
    """

    # need to check shapes of sigma_8, omega_m, likelihood
    if len(likelihood.shape)>2:
        raise NotImplementedError("Not currently prepared for likelihoods "
                                  "sampled in more than 2 dimensions")
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
    res = minimize_scalar(_compute_score, bracket = None, bounds = None,
                          args = (sigma_8, red_omega, likelihood), tol = tol,
                          options =None)
    alpha = None

    return alpha
