import unittest

import tempfile
import pandas as pd
from MWFTP.manager import run, JobInfo



class TestStringMethods(unittest.TestCase):

    def test_run(self):
        infos = [JobInfo(name='aa', info=(1, 2, 3)),
                 JobInfo(name='bb', info=(0, 0, 0))]

        def data_iterator_gen(info):
            return ((list(info), [x * x for x in info]),)

        def xy_func(x, y):
            return {'sum': sum(x) + sum(y)}

        result = []

        def collector(files):
            for f in files:
                df = pd.read_hdf(f)
                for x in df['sum'].values:
                    result.append(x)

        tmpdir = tempfile.mkdtemp()
        run(infos, data_iterator_gen, xy_func, collector, tmpdir)

        expected = [20, 0]

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
