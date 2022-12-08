import os
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_annots

# parameters
base_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/10K/within'

mwas_file_path = os.path.join(base_dir, 'mb_gwas_significant.h5')

jobs_path = os.path.join(base_dir, 'jobs')


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='20G') as q:
    q.startpermanentrun()

    snps = q.method(mwas_annots.run, (mwas_file_path, base_dir))
    q.waitforresult(snps)
