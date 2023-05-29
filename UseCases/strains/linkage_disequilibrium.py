import os
import glob
import pickle
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE

study = '10K'
min_reads_per_snp = 3

maf_file = os.path.join(f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/pcs_covariate', 'mb_snp_maf_{}.h5')
cov_file = os.path.join(f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/coverage_covariate', 'mb_snp_coverage_{}.h5')

res_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/strains/{study}/linkage_disequilibrium'


def do(species):
    maf = pd.read_hdf(maf_file.format(species))  # contains MAF_MISSING_VALUE
    contigs = set(maf.columns.str.split('_').str[2])
    for contig in contigs:
        contig_columns = maf.columns[maf.columns.str.split('_').str[2] == contig]
        contig_positions = contig_columns.str.split('_').str[0].astype(int)

        cov = pd.read_hdf(cov_file.format(species), key=f'/C_{contig}')[contig_positions]

        maf_filtered = maf.loc[cov.index, contig_columns]
        maf_filtered = maf_filtered.replace(MAF_MISSING_VALUE, np.nan)
        maf_filtered.columns = contig_positions  # necessary for mask function
        maf_filtered = maf_filtered.where(cov < min_reads_per_snp)
        maf_filtered = maf_filtered[sorted(contig_positions)]  # so pos1 will be smaller than pos2

        results = []
        for i, pos1 in enumerate(maf_filtered.columns):
            for pos2 in maf_filtered.columns[i+1:]:
                common_samples = maf_filtered[[pos1, pos2]].dropna().index
                r, p = pearsonr(maf_filtered.loc[common_samples, pos1],
                                maf_filtered.loc[common_samples, pos2])
                results.append([pos1, pos2-pos1, len(common_samples), r, p])

        with open(os.path.join(res_dir, f'{species}_C_{contig}.pkl'), 'wb') as f:
            pickle.dump(results, f)


if __name__ == '__main__':

    jobs_dir = os.path.join(os.path.dirname(res_dir), 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(res_dir)

    with qp(jobname='LD', _tryrerun=True, _delete_csh_withnoerr=True, _mem_def='20G') as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        tested = glob.glob(maf_file.format('*'))####notice this is naive tested, and not necessairly with 100 snps
        for file in tested:
            s = f'Rep_{file.split("_")[-1].split(".")[0]}'
            tkttores[s] = q.method(do, [s], _job_name=f'l{s.split("_")[-1]}')
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
