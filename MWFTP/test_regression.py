import unittest

import pandas as pd
from MWFTP.regression import regression


class TestRegression(unittest.TestCase):
    def test_regression(self):
        y = pd.DataFrame({'y': [1, 2, 3, 4]})
        x = pd.DataFrame({'a': [1, 4, 9, 16], 'b': [1, 8, 27, 64]})

        got = regression(x, y, 'my_key')

        want_keys = {
            'key', 'rsquared', 'a_coef', 'a_pval', 'a_coef_025', 'a_coef_975',
            'b_coef', 'b_pval', 'b_coef_025', 'b_coef_975'
        }
        print(got)
        self.assertEqual(got.keys(), want_keys)
        self.assertEqual(got['key'], 'my_key')


if __name__ == '__main__':
    unittest.main()
