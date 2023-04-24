import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.Utils import load_h5_files
from LabUtils.addloglevels import sethandlers

study = '10K'
run_type = 'between'  # always

run_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/{run_type}'
look_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/*_{study}_anti_mwas_{run_type}'


def rw(s, d):
    files = glob.glob(os.path.join(look_dir, 'raw_data', f'mb_gwas_Rep_*_{s}.h5'))
    if len(files) > 0:
        df = load_h5_files(files)
        df.to_hdf(os.path.join(d, 'raw_data', f'mb_gwas_Rep_all_{s}.h5'), key='snps', complevel=9)
    for file in files:
        os.remove(file)  ########notice


jobs_dir = os.path.join(run_dir, 'jobs')
os.chdir(jobs_dir)
sethandlers(file_dir=jobs_dir)

os.makedirs(os.path.join(run_dir, 'raw_data'))

with qp(jobname='concatY', _tryrerun=True) as q:
    q.startpermanentrun()
    tkttores = {}

    print('start sending jobs')
    ys = glob.glob(os.path.join(look_dir, 'raw_data', f'mb_gwas_Rep_*_Rep_*.h5'))
    ys = set(['Rep_' + file.split('_')[-1].replace('.h5', '') for file in ys])
    for species in ys:
        tkttores[species] = q.method(rw, (species, run_dir), _job_name=f'c{species.split("_")[-1]}')
    print('finished sending jobs')

    print('start waiting for jobs')
    for k, v in tkttores.items():
        q.waitforresult(v)
    print('finished waiting for jobs')

print('done')