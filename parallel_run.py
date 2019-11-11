# Python 2 (compatible with python 3)
# class to run statistics in parallel

from __future__ import absolute_import  # python2 to 3
from __future__ import division         # python2 to 3
from __future__ import print_function   # python2 to 3

try:
    import cPickle as pickle
except ImportError:
    import pickle

import os
import glob
import pandas as pd
from addloglevels import sethandlers
import shutil
import datetime


# TODO should add a more robust way to handle failures.
# Should be able to stop a run, and know which items has been finished,
# and to resume a run with a different batch size, but still run only elements that
# didn't finish yet.

class ParallelRun:
    ARG_SUF = '.pickle'
    # args_to_dump: list of args names (from func_args) that should be pickled and sent as path instead of object
    # func_args: dict of the form {arg_name: arg}, for additional arguments to "func" except for the actual elements.
    # out_path: full path to the output df
    # get_batch: generator that returns one batch of elements at a time
    def __init__(self, func, func_args,
                 out_path,
                 qp_jobname, qp_mem_def, qp_trds_def,
                 args_to_dump=None,
                 qp_kws=None, is_fake_q=False, reset_final_index=True):

        # paths
        self.out_path = os.path.abspath(out_path)
        self.run_dir = os.path.splitext(self.out_path)[0] + '.tmp'
        self._mkdirifnotexists(self.run_dir)
        self.batches_dir = os.path.splitext(self.out_path)[0] + '.batches'
        self._mkdirifnotexists(self.batches_dir)
        self.done_dir = os.path.splitext(self.out_path)[0] + '.dones'
        self._mkdirifnotexists(self.done_dir)
        self._reset_final_index = reset_final_index

        # The functionality
        self.func = func
        self.func_args = func_args.copy()
        self.args_to_dump = args_to_dump
        if self.args_to_dump is not None:
            self._dump_args()

        # qp
        self.is_fake_q = is_fake_q
        self.qp_kws = {} if qp_kws is None else qp_kws
        self.qp_kws['jobname'] = qp_jobname
        self.qp_kws['mem_def'] = qp_mem_def
        self.qp_kws['trds_def'] = qp_trds_def
        self._add_default_qp_kws()

# ------------------ INIT ------------------
    def _get_q(self):
        '''A wrapper around qp import, allowing its bypass with an argument'''
        if self.is_fake_q:
            from queue_tal.qp import fakeqp
            return fakeqp
        else:
            from queue_tal.qp import qp
            return qp

    def _add_default_qp_kws(self):
        if 'tryrerun' not in self.qp_kws:
            self.qp_kws['tryrerun'] = True
        if 'q' not in self.qp_kws:
            self.qp_kws['q'] = ['himem7.q']

    def _mkdirifnotexists(self, dirpath):
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)

    def _dump_arg(self, arg_name):
        """ Dump the object and replace it with its path in the function args
        The purpose of this functionality is to prevent sending the same argument
        to every job separately
        """
        arg_path = os.path.join(self.run_dir, arg_name + self.ARG_SUF)
        # dump the arg
        with open(arg_path, 'wb') as arg_fd:
            pickle.dump(self.func_args[arg_name], arg_fd)
        # the function argument to be the path instead of the actual object
        self.func_args[arg_name] = arg_path

    def _dump_args(self):
        if self.args_to_dump is not None:
            for arg_name in self.args_to_dump:
                self._dump_arg(arg_name)

    # ------------------ CLOSE ------------------
    def clean(self):
        self._clear_tmp_objects()
        shutil.rmtree(self.run_dir)
        shutil.rmtree(self.done_dir)
        shutil.rmtree(self.batches_dir)

    # clear the pickles that were saved
    def _clear_tmp_objects(self):
        if self.args_to_dump is not None:
            for arg_name in self.args_to_dump:
                arg_path = self.func_args[arg_name]
                os.remove(arg_path)

# ------------------ WITH compatibility ------
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # clean tmp directories only if execution completed successfully
        if exc_type is None:
            self.clean()

# ----------------- RUN ------------------
    def run(self, get_batch):
        # init SegalQueue:
        old_dir = os.getcwd()
        os.chdir(self.run_dir)
        sethandlers()

        with self._get_q()(**self.qp_kws) as q:
            q.startpermanentrun()
            self._run_in_q(q, get_batch)
        print("[ParallelRun.Run] all jobs finished, merging final df...")
        self._merge_results(get_batch.batch_size)

        # revert to normal
        os.chdir(old_dir)

    def _run_in_q(self, q, get_batch):
        waiton = []
        for batch_ind, batch in enumerate(get_batch):
            batch_path = os.path.join(self.batches_dir,
                                      'bInd{}_bSize{}.df'.format(batch_ind, get_batch.batch_size))
            if not os.path.exists(batch_path):
                waiton.append(q.method(self._apply_on_batch, (batch, batch_path)))
        q.wait(waiton)

    def _get_orig_args(self):
        orig_args = self.func_args.copy()
        for arg_name in self.args_to_dump:
            arg_path = self.func_args[arg_name]
            with open(arg_path, 'rb') as arg_fd:
                orig_args[arg_name] = pickle.load(arg_fd)
        return orig_args

    def _apply_on_batch(self, batch, batch_path):
        # batch_res = []
        batch_res = {}
        if self.args_to_dump is not None:
            orig_args = self._get_orig_args()
        else:
            orig_args = self.func_args
            
        for name, elem in batch:
            done_path = os.path.join(self.done_dir, name+'.done')
            # this element is already done:
            if os.path.exists(done_path):
                with open(done_path, 'rb') as done_f:
                    res = pickle.load(done_f)
            # should run on this element:
            else:
                res = self.func(elem, **orig_args)
                # save tmp stats to a file to mark that this element is done
                with open(done_path, 'wb') as done_f:
                    pickle.dump(res, done_f)
            # ignore items that returned None
            if res is not None:
                batch_res[name] = res
        # batch_df = pd.DataFrame(batch_res)
        batch_df = pd.DataFrame.from_dict(batch_res, orient='index')
        batch_df.to_pickle(batch_path)
        print('[ParallelRun] this batch is done {}'.format(batch_path))

    # merge the results of the batches to one df
    def _merge_results(self, batch_size):
        ### VERY slow version:
        # final_df = None
        # for batch_path in glob.glob(os.path.join(self.batches_dir,
        #                                          '*bSize{}.df'.format(batch_size))):
        #     batch_df = pd.read_pickle(batch_path)
        #     if final_df is None:
        #         final_df = batch_df
        #     else:
        #         final_df = final_df.append(batch_df)
        ## Adi: I think the following version should be much faster:
        paths = glob.glob(os.path.join(self.batches_dir, '*bSize{}.df'.format(batch_size)))
        final_df = pd.concat([pd.read_pickle(path) for path in paths])

        if self._reset_final_index:
            final_df = final_df.reset_index(drop=True)
        try:
            final_df.to_pickle(self.out_path)
            print('[ParallelRun] Saved df: {}, shape:{}'.format(self.out_path, final_df.shape))
        except: 
            tmp_path = '_'.join(str(datetime.datetime.now()).split()) + '_parallelRunTmpOut.df'
            tmp_path = os.path.join(os.path.expanduser('~'), tmp_path)
            final_df.to_pickle(tmp_path)
            print('[ParallelRun] Problem with out_path ({}), ' 
                  'Saved df to tmp path in your home folder: {}, '
                  'shape:{}'.format(self.out_path, tmp_path, final_df.shape))


class GetBatchDf:
    def __init__(self, df, batch_size, names=None):
        self.df = df
        self.batch_size = batch_size
        self.cur_row = 0
        if names is not None:
            if len(names) != self.df.shape[0]:
                raise ValueError("[GetBatchDf] names length must equal the number of rows in df")
            self.names = names
        else:
            self.names = [str(i) for i in range(self.df.shape[0])]

    def __iter__(self):
        return self

    def next(self):
        if self.cur_row >= self.df.shape[0]:
            raise StopIteration
        else:
            start = self.cur_row
            end = start + self.batch_size
            self.cur_row += self.batch_size

            cur_rows = self.df.iloc[start: end, :]
            cur_names = self.names[start:end]
            # returns an "iterable df":
            return [(name, row) for name, row in zip(cur_names, cur_rows.iterrows())]


class GetBatchIterable:
    def __init__(self, itr, batch_size, names=None):
        self.itr = itr
        self.batch_size = batch_size
        self.cur_elem = 0
        if names is not None:
            if len(names) != len(itr):
                raise ValueError("[GetBatchIterable] names and itr must have the same length")
            self.names = names
        else:
            self.names = [str(i) for i in range(len(self.itr))]

    def __iter__(self):
        return self

    def next(self):
        if self.cur_elem >= len(self.itr):
            raise StopIteration
        else:
            start = self.cur_elem
            end = start + self.batch_size
            self.cur_elem += self.batch_size

            cur_it = self.itr[start: end]
            cur_names = self.names[start: end]
            return zip(cur_names, cur_it)


if __name__ == '__main__':

    # Demo run for testing this module:
    def do_anything(num, factor, lst):
        return {'in': num,
                'out': num * factor,
                'lst': tuple(lst)}

    getter = GetBatchIterable([x for x in range(1, 11)], batch_size=3)

    another_arg = [1, 2, 3]

    out_path = r'/net/mraid08/export/jafar/Microbiome/Analyses/adiwa/metaDist/tmp/tryParallel/results.df'
    with ParallelRun(do_anything,
                     {'factor': 5, 'lst': another_arg}, out_path, 'tmpjob', '1G', 1,
                     args_to_dump=['lst'], is_fake_q=True) as runner:
        runner.run(getter)

    print("finished!")
