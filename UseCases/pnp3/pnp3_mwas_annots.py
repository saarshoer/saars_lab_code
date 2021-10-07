import os
from LabData.DataAnalyses.MBSNPs import mwas_annots
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

# parameters
body_site = 'Oral'

output_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas_{body_site.lower()}_small'
jobs_path = os.path.join(output_dir, 'jobs')
mwas_file_path = os.path.join(output_dir, f'mb_gwas_significant.h5')

# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname=f'annot_{body_site}', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=2, _mem_def='15G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annots.run, (mwas_file_path, output_dir, 'PNP3'))
    q.waitforresult(snps)
