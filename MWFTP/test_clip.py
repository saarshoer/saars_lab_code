import unittest

import pandas as pd
from MWFTP.clip import clip


class TestClip(unittest.TestCase):
    def test_clip(self):
        data = pd.Series(
            [1, -1, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        )
        got = clip(data, window_frac=1, stds_clip=1, stds_remove=2)
        want = pd.Series(
            [0.8451, -0.8451, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        )
        pd.testing.assert_series_equal(got, want, rtol=0.0001)


if __name__ == '__main__':
    unittest.main()
