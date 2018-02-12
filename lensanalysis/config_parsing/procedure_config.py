import ConfigParser

import numpy as np
from astropy import units as u

from lensanalysis.misc.name_parser import SafeNameParser


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

class ProcedureConfig(object):
    def __init__(self,procedure_config):
        assert isinstance(procedure_config,SequentialArgConfigParser)
        self._config = procedure_config

    def has_peak_count_bins(self):
        return self._config.has_section("PeakCountBins")

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
                                              "tomo_bin_min", num_bins,
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
        return cls(config_parser)
