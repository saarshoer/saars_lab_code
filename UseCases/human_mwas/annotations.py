import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_annots

# parameters
base_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/human_mwas'
maf_template = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/Cache/global_mafs/all_data/mb_snp_g_maf_{}_R1_S500.h5'

# run
jobs_path = os.path.join(base_dir, 'jobs')
os.makedirs(jobs_path, exist_ok=True)
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=False, q=['himem7.q'], _mem_def='5G') as q:
    q.startpermanentrun()

    snps = {}
    for params in ['saars_params', 'lirons_params']:
        for study in ['Lifeline_deep', '10K', 'D2']:
            for pheno in ['age', 'sex', 'bmi']:
                curr_dir = os.path.join(base_dir, params, study, pheno)
                mwas_file_path = os.path.join(curr_dir, 'mb_gwas_significant.h5')
                snps[curr_dir] = q.method(mwas_annots.run, (mwas_file_path, curr_dir, 'Gut', maf_template))

    for value in snps.values():
        q.waitforresult(value)
