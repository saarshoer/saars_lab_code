import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
from pnp3_mwas import P
output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/pull_requests/annotations/test'
mwas_file_path = '/net/mraid08/export/genie/LabData/Analyses/lironza/20200809_194907_mwas_bmi/mb_gwas.h5'
jobs_path = output_dir


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='Tannot', _delete_csh_withnoerr=False, q=['himem7.q'], max_r=1, _mem_def='10G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot.run, (P, mwas_file_path, output_dir))
    q.waitforresult(snps)
