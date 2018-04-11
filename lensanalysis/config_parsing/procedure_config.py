import ConfigParser
import os.path
from string import Formatter

import numpy as np
from astropy import units as u

from .cosmo_config import get_abs_paths
from lensanalysis.misc.name_parser import SafeNameParser
from lensanalysis.misc.fname_formatter import BaseFnameFormatter, \
    AbstractFnameFormatter


def _get_uniform_type_list(self, section,option,converter):
    temp = self.get(section,option)
    return [converter(elem) for elem in temp.split(',')]

class SequentialArgConfigParser(ConfigParser.SafeConfigParser):

    def getint_list(self,section,option):
        return _get_uniform_type_list(self, section,option,int)

    def getfloat_list(self,section,option):
        return _get_uniform_type_list(self, section,option,float)

def _build_peak_count_bins(bin_min, bin_max, bin_count,sigma):
    if bin_min>=bin_max:
        raise ValueError("bin_min must be less than bin_max")
    if bin_count <= 0:
        raise ValueError("bin_count must be at least one")
    elif bin_count == 1:
        return np.array([bin_min,bin_max]) * sigma

    return np.linspace(bin_min, bin_max, num = bin_count + 1,
                       endpoint = True) * sigma

def _load_indiv_bin_vals(config, section, option, num_bins, val_type = int):
    """
    helper function to load individual values for each tomographic bin.
    """
    if val_type == int:
        vals = config.getint_list(section,option)
    elif val_type == float:
        vals = config.getfloat_list(section,option)
    else:
        raise ValueError('val_type can only be "int" or "float".')
    if len(vals) != num_bins:
        raise ValueError(("{:d} values must be specified for the {:s} "
                          "option.").format(num_bins,option))
    return vals

def _normalize_sigma_condition_checking(config,tomo=False):
    if tomo:
        option_name = "tomo_sigma"
    else:
        option_name = "normal_sigma"
    if not config.has_option("PeakCountBins",option_name):
        return None

    if tomo:
        if config.getfloat("PeakCountBins","normal_sigma") != 0:
            raise ValueError("since 'sigma_indiv_map' is True, 'normal_sigma' "
                             "must be False")
    else:
        tomo_sigma = config.getfloat_list("PeakCountBins","tomo_sigma")
        if len(tomo_sigma) !=1 or tomo_sigma[0] != 0:
            raise ValueError("tomo_sigma must be set to a single value of 0 if "
                             "tomo_sigma_indiv_map is True")

def check_num_fields(template):
    iterable = Formatter().parse(template)

    num = 0
    for _, field_name, format_spec, _ in iterable:
        if field_name is not None:
            num +=1
    return num

class ProcedureConfig(object):
    """
    We should probably break this up into individual objects that tracks 
    configuration of different tasks.
    """

    def __init__(self,procedure_config,config_fname):
        assert isinstance(procedure_config,SequentialArgConfigParser)
        self._config = procedure_config
        self._config_fname = config_fname

    def has_peak_count_bins(self):
        return self._config.has_section("PeakCountBins")

    def sigma_indiv_map(self,tomo=False):
        """
        Whether standard deviation for SNR in peak counting is computed from 
        each individual bin. Default is False.
        """
        if tomo:
            option_name = 'tomo_sigma_indiv_map'
        else:
            option_name = 'sigma_indiv_map'
        if self._config.has_option("PeakCountBins",option_name):
            
            out = self._config.getboolean("PeakCountBins", option_name)
            if out:
                _normalize_sigma_condition_checking(self._config,tomo=tomo)
            return out
        return False

    def peak_count_bins(self, num_bins):
        if (self._config.has_option("PeakCountBins","bin_fname")
            and self._config.has_option("PeakCountBins","bin_min")):
            raise ValueError("Procedure config cannot have options for both "
                             "bin_min and bin_fname")
        elif self._config.has_option("PeakCountBins","bin_fname"):
            raise NotImplementedError("Not yet implemented")
        else:
            if self._config.has_option("PeakCountBins","normal_sigma"):
                normal_sigma = self._config.getfloat("PeakCountBins",
                                                     "normal_sigma")
                if normal_sigma == 0:
                    normal_sigma = 1.
                # implicitly checks if values are okay
                val = self.sigma_indiv_map(tomo=False)
                assert normal_sigma>0
            else:
                normal_sigma = 1.

            bin_min = self._config.getfloat("PeakCountBins","bin_min")
            bin_max = self._config.getfloat("PeakCountBins","bin_max")
            bin_count = self._config.getint("PeakCountBins","bin_count")
            return _build_peak_count_bins(bin_min, bin_max, bin_count,
                                          normal_sigma)

    def tomo_peak_count_bins(self,num_bins):
        assert isinstance(num_bins, int) and num_bins>0

        if self._config.has_option("PeakCountBins","tomo_sigma"):
            # if we wish to allow different tomographic sigmas, then we will do
            # accomplish this by specifying multiple sigma values deliminated
            # by commas
            tomo_sigma = self._config.getfloat_list("PeakCountBins",
                                                    "tomo_sigma")
            if len(tomo_sigma) not in [1,num_bins]:
                if num_bins == 1:
                    message = ("If specified, tomo_sigma must have 1 value to "
                               "be used with the 1 tomographic bin.")
                else:
                    message = ("If specified tomo_sigma must either have 1 "
                               "value to be used with all {:d}\ntomographic "
                               "bins or {:d} different "
                               "values.").format(num_bins,num_bins)
                raise ValueError(message)

            # implicitly checks that values are allowed
            sigma_indiv_map = self.sigma_indiv_map(tomo=True)
            

            for i in range(len(tomo_sigma)):
                if tomo_sigma[i] == 0:
                    tomo_sigma[i] = 1.
                assert tomo_sigma[i]>0

            if len(tomo_sigma) == 1:
                temp = tomo_sigma
                tomo_sigma = [temp[0] for i in range(num_bins)]
        else:
            tomo_sigma = [1. for i in range(num_bins)]

        const_bins = self._config.getboolean("PeakCountBins","tomo_const_bins")
        if const_bins:
            bin_max_l = []
            bin_min_l = []
            bin_count_l = []
            for i in range(num_bins):
                bin_min_l.append(self._config.getfloat("PeakCountBins",
                                                       "tomo_bin_min"))
                bin_max_l.append(self._config.getfloat("PeakCountBins",
                                                       "tomo_bin_max"))
                bin_count_l.append(self._config.getint("PeakCountBins",
                                                       "tomo_bin_count"))
        else:
            bin_min_l =_load_indiv_bin_vals(self._config, "PeakCountBins",
                                            "tomo_bin_min", num_bins,
                                            val_type = float)
            bin_max_l =_load_indiv_bin_vals(self._config, "PeakCountBins",
                                            "tomo_bin_max", num_bins,
                                            val_type = float)
            bin_count_l =_load_indiv_bin_vals(self._config, "PeakCountBins",
                                              "tomo_bin_count", num_bins,
                                              val_type = int)

        out = []
        iter_tuple = zip(bin_min_l, bin_max_l, bin_count_l, tomo_sigma)
        for bin_min, bin_max, bin_count, sigma_val in iter_tuple:
            out.append(_build_peak_count_bins(bin_min, bin_max, bin_count,
                                              sigma_val))
        return out

    def get_noise_seed(self):
        return self._config.getint("AnalysisOptions", "noise_seed")

    def diff_noise_seed_per_tomo_bin(self):

        if self._config.has_option("AnalysisOptions","separate_tomo_bin_seed"):
            return self._config.getboolean("AnalysisOptions",
                                           "separate_tomo_bin_seed")
        return True

    def noise_rs_correction(self):
        """
        Whether the rs_correction should be applied when adding noise.
        Default is False.
        """
        if self._config.has_option("AnalysisOptions","rs_correction"):
            return self._config.getboolean("AnalysisOptions", "rs_correction")
        return False

    def mask_convergence_conversion(self):
        """
        Whether the convergence map should be masked if a masked shear map is 
        converted to the convergence map. Default is False.
        """
        if self._config.has_option("AnalysisOptions",
                                   "mask_convergence_conversion"):
            return self._config.getboolean("AnalysisOptions",
                                           "mask_convergence_conversion")
        return False

    def real_space_smoothing(self):
        """
        If we should smooth the Shear Map/Convergence Map in real space (if we 
        are smoothing at all). Default is False.
        """

        if self._config.has_option("AnalysisOptions", "real_space_smoothing"):
            return self._config.getboolean("AnalysisOptions",
                                           "real_space_smoothing")
        return False

    def edge_mode(self):
        """
        Controls how smoothing handles the edge case (if we are smoothing at 
        all). Allowed values for real space smoothing are 'constant' and '
        mirror' (the latter is the default. The only option for Fourier space 
        smoothing is 'constant'.
        """
        real_space_smoothing = self.real_space_smoothing()

        if self._config.has_option("AnalysisOptions", "edge_mode"):
            val = self._config.getboolean("AnalysisOptions", "edge_mode")
            if val != 'default':
                if real_space_smoothing and val not in ['constant','mirror']:
                    raise ValueError("The only allowed values of the edge_mode "
                                     "option for real space smoothing include "
                                     "{default, constant, mirror}.")
                elif (not real_space_smoothing) and val != 'constant':
                    raise ValueError("The only allowed values of the edge_mode "
                                     "option for Fourier space smoothing "
                                     "include {default, constant}.")
                return val
        # return default values
        if real_space_smoothing:
            return "mirror"
        else:
            return "constant"
    
    def get_smoothing_scale(self):
        value = self._config.getfloat("AnalysisOptions",
                                      "smoothing_scale")
        unit = u.Unit(self._config.get("AnalysisOptions",
                                       "smoothing_unit"))
        if unit.physical_type != "angle":
            raise ValueError("smoothing_unit must be a unit in which angles "
                             "are measured")
        return value*unit

    def get_min_realization(self):
        return self._config.getint("AnalysisOptions", "min_realization")

    def get_max_realization(self):
        return self._config.getint("AnalysisOptions", "max_realization")

    def get_conv_map_resolution(self):
        return self._config.getint("ShearCatalogSpecific","conv_map_resolution")

    def get_conv_map_edge_angle(self):
        value = self._config.getfloat("ShearCatalogSpecific",
                                      "map_edge_angle")
        unit = u.Unit(self._config.get("ShearCatalogSpecific",
                                       "map_edge_angle_unit"))
        if unit.physical_type != "angle":
            raise ValueError("map_edge_angle_unit must be a unit in which "
                             "angles are measured")
        return value*unit

    def get_default_saved_products(self):

        parser = SafeNameParser()
        out = []
        items = self._config.items('AnalysisProductSaving')
        for name,value in items:
            descriptor,object_name = parser.parse_name(name)
            if object_name == "shear_cat":
                raise ValueError("shear_cat cannot be saved")
            if value:
                out.append(object_name)
        return out

    @classmethod
    def from_fname(cls,fname):
        config_parser = SequentialArgConfigParser()
        config_parser.read(fname)
        return cls(config_parser,fname)

    def has_non_ppz_rebinning(self):
        """
        Indicating whether or not there will be rebinning in the sense that the 
        bounds of one or more tomographic bins is changed. If the pseudo photoz 
        is being used, but the intervals of each tomographic remain the same 
        (and the galaxies will be rebinned so that they are binned by their 
        photoz, not their z_spec) this returns False.
        """
        result = self._config.getboolean("Rebinning","rebin")
        return result
    
    def get_bin_limits(self):
        """
        The bin limits have been explicitly specified of the resulting bins.

        If we are running analysis on a pseudo photoz dataset using the input 
        bins, this is implicitly asking for the input bins.
        """
        if self._config.has_option("Rebinning","num_bins"):
            num_bins = self._config.getint("Rebinning","num_bins")
            if num_bins == -1:
                return None
            if num_bins <1:
                raise ValueError("num_bins must be -1 or a positive integer.")
        else:
            return None

        if self._config.has_option("Rebinning","bin_limits"):
            vals = self._config.getfloat_list("Rebinning","bin_limits")
            if len(vals) == 1 and vals[0] == -1:
                return num_bins,None
            if len(vals) not in [num_bins+1, 2*num_bins]:
                raise ValueError(("bin_limits must be equal to -1 (not "
                                  "specifying any limits)\nor must provide "
                                  "either {:d} or {:d} monotonically increaing "
                                  "values (the neighboring values can also\nbe "
                                  "equal to one another.").format(num_bins+1,
                                                                  2*num_bins))

            if (np.diff(vals)<0).any():
                raise ValueError("values of bin_limits must be monotonically "
                                 "increasing (neighboring values can also\nbe "
                                 "equal to one another.")

            if len(vals) == num_bins+1:
                bin_lim = zip(vals[:-1],vals[1:])
            else:
                bin_lim = zip(vals[::2],vals[1::2])
            return num_bins,bin_lim
        return num_bins, None

    def input_bin_lim_guidelines(self):
        """
        For use with ppz without any rebinning. If individual bin_limits have 
        not been specified, returns a dictionary of guidelines for determining 
        the input bin limits.

        The keys of the dictionary include:
            "min_bin_value"
            "max_bin_value"
            "contact_intervals"
        """
        # first we make sure that the explicit bins have not been provided
        temp = self.get_bin_limits()
        if temp is not None and temp[1] is not None:
            return None

        num_specified = 0
        out = {}

        try: 
            val = self._config.getfloat("Rebinning", "min_bin_lower_bound")
            num_specified +=1
        except ConfigParser.NoOptionError:
            val = None
        out["min_bin_value"] = val


        # check to see if the upper limit should be inclusive
        try: 
            val = self._config.getfloat("Rebinning", "max_bin_upper_bound")
            num_specified +=1
        except ConfigParser.NoOptionError:
            val = None
        try:
            inclusive = self._config.getboolean("Rebinning",
                                                "max_bin_inclusive_upper_bound")
            if inclusive:
                if val is None:
                    raise ValueError("max_bin_inclusive_upper_bound should not "
                                     "be True since max_bin_upper_bound is not "
                                     "specified")
                val = np.nextafter(val,np.inf)
        except ConfigParser.NoOptionError:
            pass
        out["max_bin_value"] = val


        try: 
            val = self._config.getboolean("Rebinning", "contact_intervals")
            num_specified +=1
        except ConfigParser.NoOptionError:
            val = True
        out["contact_intervals"] = val

        if num_specified == 0:
            return None
        else:
            return out
        

    def get_z_binning_cat(self):
        """
        If applicable, this returns the directory where position Catalogs are 
        saved and the FnameFormatter for the actual file name of the position 
        catalog.

        These position catalogs are used exclusively for rebinning.
        """

        if self._config.has_option("Rebinning","z_binning_cat_fname_template"):
            if not self._config.has_option("Rebinning","z_binning_cat_dirname"):
                raise ValueError("the z_binning_cat_fname_template cannot be a "
                                 "specified option without specifying the "
                                 "z_binning_cat_dirname option")
        elif self._config.has_option("Rebinning","z_binning_cat_dirname"):
            raise ValueError("the z_binning_cat_dirname cannot be a specified "
                             "option without specifying the "
                             "z_binning_cat_fname_template option")
        else:
            return None,None

        template = self._config.get("Rebinning","z_binning_cat_fname_template")

        if template == "":
            return None,None

        # compute the number of fields in temp
        iterable = Formatter().parse(template)
        num = 0
        for _, field_name, format_spec, _ in iterable:
            if field_name is not None:
                num +=1

        fields = []
        if num>1:
            raise ValueError("The z_binning_cat_fname_template option can only "
                             "have 0 or 1 fields")
        elif num == 1:
            fields.append("bin")
        fname_formatter = BaseFnameFormatter(template,fields)
        #print fname_formatter
        #assert(isinstance(fname_formatter,AbstractFnameFormatter))

        temp = self._config.get("Rebinning","z_binning_cat_dirname")
        if temp == "":
            return os.path.dirname(self._config_fname),fname_formatter

        path = get_abs_paths(self._config,
                             [("Rebinning","z_binning_cat_dirname")],
                             os.path.abspath(self._config_fname))[0]

        return path,fname_formatter
