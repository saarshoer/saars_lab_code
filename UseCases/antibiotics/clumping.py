import os
import glob
import pickle
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

pval_col = 'Pval'
max_positions = 100
r2_cutoff = 0.5
alpha = 0.05

study = '10K'
run_type = 'between'

base_path = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/{run_type}'
sig_file = os.path.join(base_path, 'mb_gwas_significant.h5')
data_file = os.path.join(base_path, 'raw_data', 'mb_gwas_{X}_{Y}.h5' if run_type == 'within' else 'mb_gwas_Rep_all_{Y}.h5')
clump_file = os.path.join(base_path, 'clumping72', 'mb_gwas_{X}_{Y}.pkl')


def get_data(X, Y, names):

    d = pd.read_hdf(data_file.replace('{X}', X).replace('{Y}', Y))
    d = d[(d.index.get_level_values('Species') == X)]  # makes it much more memory efficient
    d = d.reset_index(names)
    d['Position'] = d['Position'].astype(int)
    d = d.groupby(names)[d.columns[-1]]

    return d


def contig_clump(contig_snps, X, Y, names):

    data = get_data(X, Y, names)

    contig_sig_snps = []
    contig_all_snps = {}

    for new_snp in contig_snps:
        add = True
        for sig_snp in contig_sig_snps:
            if (new_snp[2].split('_')[1] == sig_snp[2].split('_')[1]) & (
                    abs((new_snp[3] - sig_snp[3])) <= max_positions):
                common_samples = set(data.get_group(new_snp).index) & set(data.get_group(sig_snp).index)
                r, p = pearsonr(data.get_group(new_snp).loc[common_samples],
                                data.get_group(sig_snp).loc[common_samples])
                if (r * r >= r2_cutoff) & (p < alpha):
                    add = False
                    break
        if add:
            contig_sig_snps.append(new_snp)
            contig_all_snps[new_snp] = new_snp
        else:
            contig_all_snps[new_snp] = sig_snp

    return contig_sig_snps, contig_all_snps


def clump(X, Y):

    sig_df = pd.read_hdf(sig_file)
    sig_df = sig_df[(sig_df.index.get_level_values('Species') == X) & (sig_df.index.get_level_values('Y') == Y)]
    # sig_df = pd.concat([sig_df.head(100), sig_df.tail(5)])
    sig_df = sig_df.sort_values(pval_col, ascending=True)

    sig_snps = []
    all_snps = {}

    with qp(jobname='RBKNcSlumping', _tryrerun=True, _mem_def='20G') as q:
        q.startpermanentrun()
        tkttores = {}

        for contig, contig_df in sig_df.groupby('Contig', as_index=False):
            tkttores[contig] = q.method(contig_clump, (contig_df.index, X, Y, sig_df.index.names))

        for k, v in tkttores.items():
            contig_sig_snps, contig_all_snps = q.waitforresult(v)

            sig_snps = sig_snps + contig_sig_snps
            all_snps.update(contig_all_snps)

    sig_snps = sig_df.loc[sig_snps].sort_values(pval_col, ascending=True).index.tolist()

    data = get_data(X, Y, sig_df.index.names)

    for i, snp2 in enumerate(sig_snps[::-1]):  # from weak to strong
        for snp1 in sig_snps[:-i-1]:  # from strong to weak
            common_samples = set(data.get_group(snp1).index) & set(data.get_group(snp2).index)
            r, p = pearsonr(data.get_group(snp1).loc[common_samples], data.get_group(snp2).loc[common_samples])
            if (r * r >= r2_cutoff) & (p < alpha):
                for snp in all_snps.keys():
                    if all_snps[snp] == snp2:
                        all_snps[snp] = snp1

    results_file = clump_file.replace('{X}', X).replace('{Y}', Y)
    with open(results_file, 'wb') as f:
        pickle.dump(all_snps, f)

    return results_file


if __name__ == '__main__':

    # queue
    jobs_dir = os.path.join(base_path, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(os.path.dirname(clump_file))#, exist_ok=True)

    with qp(jobname='BKclumping', _tryrerun=True, _mem_def='100G') as q:
        q.startpermanentrun()
        tkttores = {}

        significant_df = pd.read_hdf(sig_file)
        # significant_df = significant_df[significant_df.index.get_level_values('Species') == 'Rep_595']

        print('start sending jobs')
        runs = significant_df.reset_index()[['Species', 'Y']].drop_duplicates().reset_index(drop=True)
        for i, (speciesX, speciesY) in runs.iterrows():
            if not os.path.exists(clump_file.replace('72', '7').replace('{X}', speciesX).replace('{Y}', speciesY)):
                tkttores[i] = q.method(clump, (speciesX, speciesY))
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

        print('start df update')
        results_files = glob.glob(clump_file.replace('{X}', '*').replace('{Y}', '*'))
        significant_df['clumping'] = None
        for results_file in results_files:
            with open(results_file, 'rb') as f:
                run_results = pickle.load(f)
            significant_df.loc[run_results.keys()] = significant_df.loc[run_results.keys()].assign(
                clumping=run_results.values())
        significant_df.to_hdf(sig_file.replace('.h5', '_clumping72.h5'), key='snps')
        print('finished df update')
