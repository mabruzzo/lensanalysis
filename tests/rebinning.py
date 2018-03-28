from copy import deepcopy
import unittest

import numpy as np
from lenstools.catalog import ShearCatalog

from lensanalysis.procedure.tomo_binning import PseudoPhotozRebinner, \
    ConstantBias, build_static_z_rebinner

def _build_simple_catalogs(z_arrays,seed=None):
    np.random.seed(seed)
    out = []
    for elem in z_arrays:
        z_array = np.array(elem)
        length = np.alen(z_array)
        shear1 = np.random.uniform(0.0,1.0,(length,))
        shear2 = np.random.uniform(0.0,1.0,(length,))

        x = np.random.uniform(0.0,3.5,(length,))
        y = np.random.uniform(0.0,3.5,(length,))

        cat = ShearCatalog([shear1,shear2,z_array,x,y],
                           names = ("shear1","shear2","z","x","y"))
        out.append(cat)
    return out

def _compare_nonz_catalog_rows(cat1,cat2,rows1,rows2):
    colnames = ["shear1","shear2","x","y"]
    w_rows1 = (np.array(rows1),)
    w_rows2 = (np.array(rows2),)

    return (cat1[colnames][w_rows1] == cat2[colnames][w_rows2]).all()

class SpecificBinsPPZRebinnerConstantBiasTestCase(unittest.TestCase):
    def setUp(self):
        noise_function = ConstantBias(0.003)
        self.rebinner = PseudoPhotozRebinner(noise_function,[(0,1),(1,2)])
        #self.neg_rebinner = PseudoPhotozRebinner(ConstantBias(-0.003),
        #                                         [(0,1),(1,2)])
        self.larger_rebinner = PseudoPhotozRebinner(noise_function,
                                                    [(0,1),(1.003,2),(2.2,3)])

    def test_basic_catalog(self):
        """
        This is a very straight-forward test.
        """
        temp = _build_simple_catalogs([[0.00, 0.5,0.90, 0.992, 0.995, 0.999],
                                       [1.5,1.6]],
                                      seed=57)
        #print temp
        data_objects = deepcopy(temp)
        result = self.rebinner.rebin(data_objects,3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[0],[0,1,2,3],
                                                   [0,1,2,3]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[1],[4,5],
                                                   [0,1]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[1],result[1],
                                                   [0,1],[2,3]))

    def test_overflow_basic_catalog(self):
        """
        This test test checks that a galaxy with a photo-z outside of the 
        range correctly get's clipped from the resulting catalogs.
        """
        temp = _build_simple_catalogs([[0.00, 0.5,0.90, 0.992, 0.995, 0.999],
                                       [1.5,1.6,1.992]],
                                      seed=57)
        #print temp
        data_objects = deepcopy(temp)
        result = self.rebinner.rebin(data_objects,3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[0],[0,1,2,3],
                                                   [0,1,2,3]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[1],[4,5],
                                                   [0,1]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[1],result[1],
                                                   [0,1],[2,3]))
        self.assertTrue(len(result[1]) == 4)

    def test_between_intervals_more_bins(self):
        """
        Using more intervals, we test to see that the photometrically shifted 
        redshifts lying between intervals are excluded from the tomographic 
        bins.
        """
        temp = _build_simple_catalogs([[0.00, 0.5,0.90, 0.992, 0.995, 0.999],
                                       [1.5,1.6,1.992],
                                       [2.7,2.8]],
                                      seed=57)
        #print temp
        data_objects = deepcopy(temp)
        result = self.larger_rebinner.rebin(data_objects,3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[0],[0,1,2,3],
                                                   [0,1,2,3]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[1],[5],
                                                   [0]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[1],result[1],
                                                   [0,1],[1,2]))
        self.assertTrue(len(result[1]) == 3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[2],result[2],
                                                  [0,1],[0,1]))
        self.assertTrue(len(result[2]) == 2)

    def test_overflow_more_bins(self):
        """
        we are testing the case where bias pushes something from the middle 
        bin beyond the upper bin.
        """
        rebinner = PseudoPhotozRebinner(ConstantBias(0.003),
                                        [(0,1),(1.003,2),(2,2.003)])

        temp = _build_simple_catalogs([[0.00, 0.5,0.90, 0.992, 0.995, 0.999],
                                       [1.5,1.6,1.992,1.995,1.999],
                                       [2.02]],
                                      seed=57)

        data_objects = deepcopy(temp)
        result = rebinner.rebin(data_objects,3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[0],[0,1,2,3],
                                                   [0,1,2,3]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[0],result[1],[5],
                                                   [0]))
        self.assertTrue(_compare_nonz_catalog_rows(temp[1],result[1],
                                                   [0,1],[1,2]))
        self.assertTrue(len(result[1]) == 3)
        self.assertTrue(_compare_nonz_catalog_rows(temp[1],result[2],
                                                   [2],[0]))
        self.assertTrue(len(result[2]) == 1)
        
class StaticRebinnerTestCase(unittest.TestCase):
    def test_subdivide_rebin(self):
        # going to start with 2 catalogs and wind up with 4 catalogs
        z_arrays = [[0.25,0.356,0.536,0.8536],
                    [1.163,1.325,1.634,1.934]]
        rebin_z_intervals = [(0.,0.52),(0.52,1.1),(1.1,1.45),(1.45,2.0)]

        data_object = _build_simple_catalogs(z_arrays,seed=55)
        rebinner = build_static_z_rebinner(data_object, rebin_z_intervals,
                                           share_position_component = False)

        rebinned_do = rebinner.rebin(data_object,0)
        #cat2 = _build_simple_catalogs(z_arrays,seed=55)

        for i in [0,1]:
            # first 2 rows of original ith cat should be both rows in the
            # rebinned (2*i)th cat 
            self.assertTrue(_compare_nonz_catalog_rows(data_object[i],
                                                       rebinned_do[2*i],
                                                       [0,1],[0,1]))
            self.assertTrue(len(rebinned_do[2*i]) == 2)
            # second 2 rows of original ith cat should be both rows in the
            # rebinned (2*i+1)th cat
            self.assertTrue(_compare_nonz_catalog_rows(data_object[i],
                                                       rebinned_do[2*i+1],
                                                       [2,3],[0,1]))
            self.assertTrue(len(rebinned_do[2*i+1]) == 2)

    def test_subdivide_rebin_shared_col(self):
        # going to start with 2 catalogs and wind up with 4 catalogs
        z_arrays = [[0.25,0.356,0.536,0.8536],
                    [1.163,1.325,1.634,1.934]]
        rebin_z_intervals = [(0.,0.52),(0.52,1.1),(1.1,1.45),(1.45,2.0)]

        data_object = _build_simple_catalogs(z_arrays,seed=55)
        rebinner = build_static_z_rebinner(data_object, rebin_z_intervals,
                                           share_position_component = True)

        rebinned_do = rebinner.rebin(data_object,0)
        data_object2 = _build_simple_catalogs(z_arrays,seed=59)

        pass
        

def suite():
    tests = ['test_basic_catalog','test_overflow_basic_catalog',
             'test_between_intervals_more_bins',
             'test_overflow_more_bins']
    test_case_class = SpecificBinsPPZRebinnerTestCase
    s = unittest.TestSuite()
    s.addTests(map(test_case_class, tests))

    tests = ['test_subdivide_rebin','test_subdivide_rebin_shared_col']
    test_case_class = StaticRebinnerTestCase
    s.addTests(map(test_case_class, tests))
    return s

if __name__ == '__main__':
    unittest.main()
