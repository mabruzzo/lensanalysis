import numpy as np
from lenstools import ShearMap, ConvergenceMap

from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

def _determine_mask(shear_map,map_mask):
    assert map_mask is None or isinstance(map_mask,np.ndarray)
    if map_mask is None:
        temp = np.isnan(shear_map.data)
        mask = np.logical_or(temp[0,:,:], temp[1,:,:])
    else:
        assert map_mask.shape == shear.data.shape[1:]
        mask = map_mask
    return mask

def _unmasking(shear_map,mask,fill):
    assert isinstance(fill,(int,float))
    shear_map.data[0][mask] = 0.0
    shear_map.data[1][mask] = 0.0

def _remasking(shear_map,conv_map,mask,mask_conv = True):
    if mask_conv:
        conv_map.mask(np.logical_not(mask).astype(np.int8),
                      inplace = True)
    shear_map.data[0][mask] = np.nan
    shear_map.data[1][mask] = np.nan

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

    mask = _determine_mask(shear_map,map_mask)
    _handle_unmasking(shear_map,mask,fill)
    conv_map = shear_map.convergence()
    _remasking(shear_map,conv_map,mask)
    return conv_map

def _get_smooth_scale(map,scale_angle):
    assert (scale_angle.unit.physical_type ==
            map.side_angle.unit.physical_type)
    # convert smoothing_scale from units of angle to units of pixels.
    sigma_pix = (scale_angle * map.data.shape[-1] 
                 / map.side_angle).decompose().value
    assert sigma_pix > 0
    return sigma_pix

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

    sigma_pix = _get_smooth_scale(conv_map,scale_angle)
    
    mask = np.isnan(conv_map.data)

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
        if hasattr(conv_map,"gradient_x") or hasattr(conv_map,"gradient_y"):
	    conv_map.gradient()

	if (hasattr(conv_map,"hessian_xx") or hasattr(conv_map,"hessian_yy")
            or hasattr(conv_map,"hessian_xy")):
	    conv_map.hessian()
			
        return None
    else:
        kwargs = dict()
        for attribute in conv_map._extra_attributes:
            kwargs[attribute] = getattr(conv_map,attribute)

        temp =  conv_map.__class__(temp, conv_map.side_angle,
                                   masked = preserve_mask, **kwargs)
        if preserve_mask:
            temp.mask(np.logical_not(mask).astype(np.int8),
                      inplace = True)
        return temp


def _zero_padded_fourierEB(shear_map, padding):
    data = np.lib.pad(shear_map.data,[(0,0),(0,padding),(0,padding)],
                      'constant')
    # This following is taken from the convergence instance method of the
    # ShearMap class defined in lenstools

    ft_data1 = np.fft.rfft2(data[0])
    ft_data2 = np.fft.rfft2(data[1])

    #Compute frequencies
    lx = np.fft.rfftfreq(ft_data1.shape[0])
    ly = np.fft.fftfreq(ft_data1.shape[0])

    #Safety check
    assert len(lx)==ft_data1.shape[1]
    assert len(ly)==ft_data1.shape[0]

    #Compute sines and cosines of rotation angles
    l_squared = lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2
    l_squared[0,0] = 1.0

    sin_2_phi = 2.0 * lx[np.newaxis,:] * ly[:,np.newaxis] / l_squared
    cos_2_phi = (lx[np.newaxis,:]**2 - ly[:,np.newaxis]**2) / l_squared

    #Compute E and B components
    ft_E = cos_2_phi * ft_data1 + sin_2_phi * ft_data2
    ft_B = -1.0 * sin_2_phi * ft_data1 + cos_2_phi * ft_data2

    ft_E[0,0] = 0.0
    ft_B[0,0] = 0.0

    assert ft_E.shape == ft_B.shape
    assert ft_E.shape == ft_data1.shape

    return ft_E,ft_B

def convert_shear_to_smoothed_convergence_main(shear_map, pad_axis,
                                               fft_kernel = None,
                                               map_mask = None, fill = 0,
                                               mask_result = False):
    """
    Main function for converting shear map to smoothed convergence map.

    This function assumes that the smoothing kernel has already been 
    transformed.

    Parameters
    ----------
    shear_map : lenstools.ShearMap
        The shear map we will be converting. If mask is None, then I will 
        just assume everywhere that is NaN, is to be masked.
    pad_axis : int
        The number pixels that the shear map be zero-padded along axes 1 and 2 
        (to avoid circular convolution.)
    fft_kernel : np.ndarray, None
        The fouier transformed smoothing kernel. If None, then the convergence 
        map is not smoothed. Default is None.
    mask : np.ndarray,optional
        Default is None. If this is not None, then we assume all values of NaN 
        in the shear map are masked. This must be the same size as shear_map
    fill : float or int, optional
        The value we fill in at the masked value when we perform the fourier 
        transform. Default is 0.
    mask_result : bool, optional
        Whether or not the resulting convergence map should be masked. Default 
        is False.
    """
    assert fft_kernel is None or isinstance(fft_kernel,np.ndarray)
    assert pad_axis >= 0
    # first we identify the mask and unmask the shear map
    mask = _determine_mask(shear_map,map_mask)
    _unmasking(shear_map,mask,fill)

    # now get the fourier transformed E and B mode of the shear map
    ft_E,ft_B = _zero_padded_fourierEB(shear_map, pad_axis)

    if fft_kernel is None:
        final_E = ft_E
    else:
        final_E = ft_E * fft_kernel

    # transform back
    temp = np.fft.irfft2(final_E)
    axis0,axis1 = temp.shape
    kwargs = dict((k,getattr(self,k)) for k in shear_map._extra_attributes)
    conv_map = ConvergenceMap(temp[:(axis0-pad_axis),:(axis1-pad_axis)],
                              shear_map.side_angle, **kwargs)
    _remasking(shear_map,conv_map,mask,mask_conv =mask_result)
    return conv_map


def determine_kernel_fft(sigma_pix, conv_map_length, truncate = 4.0):
    """
    Computes the fft of gaussian kernel and computes amount of padding.

    Parameters
    ----------
    sigma_pix : float,int
        The standard deviation of the Gaussian Kernel in units of pixels
    conv_map_length : int
        The number of pixels along each axis of the unpadded convergence map.
    truncate : float
        Truncate the Gaussian Filter at this many standard deviations
    """
    if truncate == 4.0:
        kernel = Gaussian2DKernel(sigma_pix).array
        kernel/=np.sum(kernel)
    else:
        raise NotImplementedError("Not currently equipped to handle other "
                                  "truncation sizes")

    pad_axis = kernel.shape[0]
    new_shape = np.array([pad_axis+conv_map_length for i in range(2)])
    if new_shape[0]%2 != 0:
        # things get weird when the shape is not even
        new_shape[0] +=1
        new_shape[1] +=1
        pad_axis+=1
    kern_shape = kernel.shape

    # taken from astropy's convolve fft
    kernslices = []
    for newdimsize, kerndimsize in zip(new_shape,kern_shape):
        center = newdimsize - (newdimsize + 1) // 2
        kernslices.append(slice(center - kerndimsize // 2,
                                center + (kerndimsize + 1) // 2))
    big_kernel = np.zeros(new_shape,dtype=kernel.dtype)
    big_kernel[kernslices] = kernel
    # resize the kernel - the kernel is now centered in big_kernel
    # now, shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    shifted_kernel = np.fft.ifftshift(big_kernel)

    # need to to take fft of kernel
    return np.fft.rfft2(shifted_kernel),pad_axis


def convert_shear_to_smoothed_convergence(shear_map, scale_angle,
                                          truncate = 4.0, map_mask = None,
                                          fill = 0, mask_result = False):
    """
    Converting shear map to smoothed convergence map.

    Parameters
    ----------
    shear_map : lenstools.ShearMap
        The shear map we will be converting. If mask is None, then I will 
        just assume everywhere that is NaN, is to be masked.
    scale_angle : astropy.Quantity
        The anglular size of one standard deviation of the Gaussian.
    truncate : float
        Truncate the Gaussian Filter at this many standard deviations
    mask : np.ndarray,optional
        Default is None. If this is not None, then we assume all values of NaN 
        in the shear map are masked. This must be the same size as shear_map
    fill : float or int, optional
        The value we fill in at the masked value when we perform the fourier 
        transform. Default is 0.
    mask_result : bool, optional
        Whether or not the resulting convergence map should be masked. Default 
        is False.
    """
    sigma_pix = _get_smooth_scale(shear_map,scale_angle)
    conv_map_length = shear_map.data.shape[-1]
    fft_kernel,pad_axis = determine_kernel_fft(sigma_pix, conv_map_length,
                                               truncate)

    return convert_shear_to_smoothed_convergence_main(shear_map, pad_axis,
                                                      fft_kernel, map_mask,
                                                      fill, mask_result)
