import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.Utils import load_h5_files
from LabUtils.addloglevels import sethandlers

study = 'Lifeline_deep'
# run_type = 'saars_params'
run_type = 'lirons_params'

run_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/human_mwas/{run_type}/{study}'
look_dir = f'/net/mraid20/export/genie/LabData/Analyses/saarsh/20250207_1829*_{study}_'


def rw(s, d):
    files = glob.glob(os.path.join(look_dir+s, 'raw_data', f'mb_gwas_Rep_*_{s}.h5'))
    if len(files) > 0:
        df = load_h5_files(files)
        df.to_hdf(os.path.join(d, s, f'mb_gwas_data.h5'), key='snps', complevel=9)


jobs_dir = os.path.join(os.path.dirname(os.path.dirname(run_dir)), 'jobs')
os.chdir(jobs_dir)
sethandlers(file_dir=jobs_dir)

with fakeqp(jobname='concatY', _tryrerun=True) as q:
    q.startpermanentrun()
    tkttores = {}

    print('start sending jobs')
    ys = glob.glob(os.path.join(look_dir+'*', 'raw_data', f'mb_gwas_Rep_*_*.h5'))
    ys = set([file.split('_')[-1].replace('.h5', '') for file in ys])
    for species in ys:
        tkttores[species] = q.method(rw, (species, run_dir), _job_name=f'c{species.split("_")[-1]}')
    print('finished sending jobs')

    print('start waiting for jobs')
    for k, v in tkttores.items():
        q.waitforresult(v)
    print('finished waiting for jobs')

print('done')
