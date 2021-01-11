import os
import mwas_annot2
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
from pnp3_mwas import P
output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/all'
mwas_file_path = os.path.join(output_dir, 'significant_snps.h5')
jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/jobs'


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='5G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot2.run, (P, mwas_file_path, output_dir))
    q.waitforresult(snps)
