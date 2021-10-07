import os
import mwas_annots
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# parameters
output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/annotations'

x_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_raw/mb_gwas_SGB_*.h5'
y_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/*/SGB_*.h5'

mwas_file_path = os.path.join(output_dir, 'snps_codons.h5')

jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs'


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='20G') as q:
    q.startpermanentrun()

    # snps_unique = q.method(mwas_annots.find_unique_snps,
    #                        (x_mwas_files_path, y_mwas_files_path, output_dir, 'Pval', 0.05/26068850133))
    # q.waitforresult(snps_unique)

    snps = q.method(mwas_annots.run, (mwas_file_path, output_dir))
    q.waitforresult(snps)
