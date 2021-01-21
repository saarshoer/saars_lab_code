import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/all'
mwas_file_path = os.path.join(output_dir, 'significant_snps_oral.h5')
jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/jobs'


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=2, _mem_def='5G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot.run, (mwas_file_path, output_dir, 'Oral'))
    q.waitforresult(snps)
