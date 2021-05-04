import os
import pandas as pd
from LabQueue.qp import qp, fakeqp
from typing import Iterable, Callable, Dict, NamedTuple, Any
from LabUtils.addloglevels import sethandlers


class JobInfo(NamedTuple):
    name: str = ''
    info: Any = None

    pass


def run_per_job(job_info: JobInfo, data_iterator_gen: Callable, xy_function: Callable[..., Dict], output_dir: str):
    results = []

    for x, y in data_iterator_gen(job_info.info):
        results.append(xy_function(x, y))

    result_df = pd.DataFrame(results)
    output_fname = os.path.join(output_dir, f'{job_info.name}.h5')
    # TODO: Investigate what the best output format is.
    result_df.to_hdf(output_fname, key='default')

    return output_fname


def run(job_iterator: Iterable[JobInfo], data_iterator: Callable, xy_function: Callable, collector: Callable,
        output_dir: str, qp_kwargs: Dict = None):
    # os.chdir('.')
    sethandlers()
    if qp_kwargs is None:
        qp_kwargs = {}

    with fakeqp(jobname='manager', **qp_kwargs) as q:
        q.startpermanentrun()

        tkttores = []

        for job_info in job_iterator:
            tkttores.append(q.method(run_per_job, (job_info, data_iterator, xy_function, output_dir),
                                    _job_name=job_info.name))

        results = []
        for r in tkttores:
            results.append(q.waitforresult(r))

        collector(results)


if __name__ == '__main__':
    qp_kwargs = {'_delete_csh_withnoerr': True, '_tryrerun': True}
    run(qp_kwargs)
