import numpy as np
from lenstools import ShearMap

from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

def convert_shear_to_convergence(shear_map, map_mask=None, fill = 0):
    """
    This function defines how we will convert from the shear map to the 
    convergence map.

    Parameters
    ----------
    shear_map : lenstools.ShearMap
        The shear map we will be converting. If mask is None, then I will 
        just assume everywhere that is NaN, is to be masked.
    mask : np.ndarray,optional
        Default is None. If this is not None, then we assume all values of NaN 
        in the shear map are masked. This must be the same size as shear_map
    fill : float or int, optional
        The value we fill in at the masked value when we perform the fourier 
        transform. Default is 0.
    """

    assert map_mask is None or isinstance(map_mask,np.ndarray)
    if map_mask is None:
        temp = np.isnan(shear_map.data)
        mask = np.logical_or(temp[0,:,:], temp[1,:,:])
    else:
        assert map_mask.shape == shear.data.shape[1:]
        mask = map_mask

    assert isinstance(fill,(int,float))
    shear_map.data[mask] = 0.0

    conv_map = shear_map.convergence()
    conv_map.mask(np.logical_not(mask).astype(np.int8),
                  inplace = True)
    shear_map[mask] = np.nan
    return conv_map

def smooth_conv_map(conv_map,scale_angle,kind="gaussian", truncate = 4.0,
                    nan_treatment = 'interpolate', inplace = False,
                    preserve_mask = True, **kwargs):
    """
    Designed for smoothing convergence maps that are masked.

    Parameters
    ----------
    conv_map : lenstools.ConvergenceMap
        The convergence map to be smoothed
    scale_angle : astropy.Quantity
        The angle of the Gaussian
    kind : str
        Type of smoothing to be performed. Select "gaussian" for regular 
        Gaussian smoothing in real space or "gaussianFFT" if you want the 
        smoothing to be performed via FFTs (advised for large scale_angle)
    truncate : float
        Truncate the Gaussian Filter at this many standard deviations
    nan_treatment : 'interpolate', 'fill'
        Whether to fill pixels or interpolate the kernel
    inplace : bool
        If True, performs the smoothing in place overwriting the old 
        convergence map
    preserve_mask : bool
    """


    assert scale_angle.unit.physical_type==self.side_angle.unit.physical_type

    mask = np.isnan(conv_map.data)
    # convert smoothing_scale from units of angle to units of pixels.
    sigma_pix = (scale_angle
                 * (conv_map.data.shape[0] / conv_map.side_angle)).value
    assert sigma_pix > 0

    if truncate == 4.0:
            kernel = Gaussian2DKernel(sigma_pix)
        else:
            raise NotImplementedError("Not currently equipped to handle other "
                                      "truncation sizes")
    if kind == 'gaussian':
        temp = convolve(conv_map.data,kernel,nan_treatment = nan_treatment,
                        **kwargs)

    elif kind == 'gaussianfft':
        temp = convolve_fft(conv_map.data,kernel,nan_treatment = nan_treatment,
                            **kwargs)

    if not conv_map._masked:
        preserve_mask = False

    if preserve_mask:
        temp[mask] = np.nan
        
    if inplace:
        conv_map.data = temp
        if hasattr(conv_map,"gradient_x") or hasattr(conv_map,"gradient_y")):
	    conv_map.gradient()

	if ((hasattr(conv_map,"hessian_xx") or hasattr(conv_map,"hessian_yy")
             or hasattr(conv_map,"hessian_xy")):
	    conv_map.hessian()
			
        return None
    else:
        kwargs = dict()
        for attribute in conv_map._extra_attributes:
            kwargs[attribute] = getattr(conv_map,attribute)

        return conv_map.__class__(temp, conv_map.side_angle,
                                  masked = preserve_mask, **kwargs)

