import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_annots

# parameters
base_dir = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/10K/within'
maf_template = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/Cache/global_mafs/all_data/mb_snp_g_maf_{}_R1_S500.h5'

mwas_file_path = os.path.join(base_dir, 'mb_gwas_significant.h5')
# mwas_file_path = os.path.join(base_dir, 'raw_hdfs', 'mb_gwas_Rep_*_Rep_*.h5')

jobs_path = os.path.join(base_dir, 'jobs')


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=False, q=['himem7.q'], _mem_def='20G') as q:
    q.startpermanentrun()

    snps = q.method(mwas_annots.run, (mwas_file_path, base_dir, 'Gut', maf_template))
    q.waitforresult(snps)

    # snps = {}
    # all_files = glob.glob(mwas_file_path)
    # for file in all_files:
    #     snps[file] = q.method(mwas_annots.run, (file, os.path.join(base_dir, 'full_annotations'), 'Gut', maf_template))
    #
    # for value in snps.values():
    #     q.waitforresult(value)
