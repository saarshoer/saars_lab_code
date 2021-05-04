import unittest

import tempfile
import pandas as pd
from MWFTP.manager import run, JobInfo


class TestManager(unittest.TestCase):

    def test_run(self):
        infos = [JobInfo(name='aa', info=(1, 2, 3)),
                 JobInfo(name='bb', info=(0, 0, 0))]

        def data_iterator_gen(info):
            return (list(info), [x * x for x in info]),

        def xy_func(x, y):
            return {'sum': sum(x) + sum(y)}

        tmpdir = tempfile.mkdtemp()
        result = run(infos, data_iterator_gen, xy_func, tmpdir)
        expected = pd.DataFrame({'sum': [20, 0]})

        pd.testing.assert_frame_equal(result, expected)


if __name__ == '__main__':
    unittest.main()
