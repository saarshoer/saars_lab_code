import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

# parameters
body_site = 'Gut'  # TODO if Oral don't forget to update majmin

output_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/PNP3_mwas_{body_site.lower()}_0months_log2division'
jobs_path = os.path.join(output_dir, 'jobs')
mwas_file_path = os.path.join(output_dir, f'mb_gwas_significant.h5')

# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname=f'annot_{body_site}', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=2, _mem_def='5G') as q:
    q.startpermanentrun()
    snps = q.method(mwas_annot.run, (mwas_file_path, output_dir, body_site))
    q.waitforresult(snps)
