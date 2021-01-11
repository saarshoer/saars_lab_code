import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
from anti_mwas import P
output_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/all/PNP3_gut_mwas_0months_delta_species_exist'

x_mwas_files_path = os.path.join(output_path, 'mb_gwas.h5')

jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/jobs'


# run
os.chdir(jobs_path)
sethandlers()

with fakeqp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='20G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot.run, (P, output_path, x_mwas_files_path))
    q.waitforresult(snps)
