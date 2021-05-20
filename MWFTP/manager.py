"""Main entry point for running the MWFTP module."""
import os
from typing import Iterable, Callable, Dict, NamedTuple, Any

import pandas as pd
from LabQueue.qp import fakeqp, qp
from LabUtils.addloglevels import sethandlers


class JobInfo(NamedTuple):
    """The information passed to a single job, to indicate its name and what data it should operate
    on."""
    name: str = ''
    info: Any = None


def _run_per_job(job_info: JobInfo, data_iterator_gen: Callable, xy_function: Callable[..., Dict],
                 output_dir: str):
    """The function that is called by each job."""
    results = []

    for x, y, k in data_iterator_gen(job_info.info):
        results.append(xy_function(x, y, k))

    result_df = pd.DataFrame(results)
    output_fname = os.path.join(output_dir, f'{job_info.name}.h5')
    # TODO: Investigate what the best output format is.
    result_df.to_hdf(output_fname, key='default')

    return output_fname


def run(job_iterator: Iterable[JobInfo], data_iterator: Callable,
        xy_function: Callable, output_dir: str, use_fakeqp=False,
        qp_kwargs: Dict = None) -> pd.DataFrame:
    """Creates a job for each item in job_iterator and collects the results."""
    sethandlers()
    if qp_kwargs is None:
        qp_kwargs = {}
    qp_defaults = {'jobname': 'manager'}

    qprovider = qp if not use_fakeqp else fakeqp

    with qprovider(**{**qp_defaults, **qp_kwargs}) as q:
        q.startpermanentrun()

        tkttores = []
        for job_info in job_iterator:
            tkttores.append(q.method(_run_per_job, (job_info, data_iterator, xy_function,
                                                    output_dir)))

        fnames = []
        for r in tkttores:
            fnames.append(q.waitforresult(r))

    result = pd.concat((pd.read_hdf(f) for f in fnames), ignore_index=True)
    for f in fnames:
        os.remove(f)

    return result
