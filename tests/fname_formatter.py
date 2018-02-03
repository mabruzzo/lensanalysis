import unittest
from lensanalysis.misc.fname_formatter import \
    RealizationBinnedSubdirectoryFormatter, BaseFnameFormatter

class SimpleRealizationBinnedSubdirectoryFormatterTestCase(unittest.TestCase):

    def setUp(self):

        # should probably mock BaseFnameFormatter
        temp = BaseFnameFormatter("WLshear_cat_bin{:d}_{:04d}r.fits",
                                  ("bin_num","realization"))
        temp = RealizationBinnedSubdirectoryFormatter('{:d}-{:d}',
                                                      ('min','max'),
                                                      288,
                                                      temp,
                                                      'realization')
        self.formatter = temp

    def test_fields(self):
        self.assertEqual(self.formatter.fields,
                         ['bin_num', 'realization'])

    def test_simple_input(self):
        self.assertEqual(self.formatter.format_fname(**{"realization":27,
                                                        "bin_num" : 4}),
                         '1-288/WLshear_cat_bin4_0027r.fits')

    def test_max_realization_in_bin(self):
        self.assertEqual(self.formatter.format_fname(**{"realization":288,
                                                        "bin_num" : 2}),
                         '1-288/WLshear_cat_bin2_0288r.fits')

    def test_min_realization_second_bin(self):
        self.assertEqual(self.formatter.format_fname(**{"realization":289,
                                                        "bin_num" : 5}),
                         '289-576/WLshear_cat_bin5_0289r.fits')

    def test_determine_subdir_1(self):
        self.assertEqual(self.formatter.determine_subdir(1,1000),
                         [('1-288',None),
                          ('289-576',None),
                          ('577-864',None),
                          ('865-1152',None)])

    def test_determine_subdir_2(self):
        self.assertEqual(self.formatter.determine_subdir(589,589),
                         [('577-864',None)])

    def test_determine_subdir_3(self):
        self.assertEqual(self.formatter.determine_subdir(352,900),
                         [('289-576',None),
                          ('577-864',None),
                          ('865-1152',None)])

def suite():
    tests = ['test_simple_input', 'test_max_realization_in_bin',
             'test_min_realization_second_bin', 'test_determine_subdir_1',
             'test_determine_subdir_2', 'test_determine_subdir_3']
    test_case_class = SimpleRealizationBinnedSubdirectoryFormatterTestCase
    return unittest.TestSuite(map(test_case_class, tests))

if __name__ == '__main__':
    unittest.main()
