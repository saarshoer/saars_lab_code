import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
from anti_mwas import P
output_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed'

x_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_raw/mb_gwas_SGB_*.h5'
y_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/*/SGB_*.h5'

jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs'


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='50G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot.run, (P, output_path, x_mwas_files_path, y_mwas_files_path))
    q.waitforresult(snps)
