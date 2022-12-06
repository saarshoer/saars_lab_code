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
abs_corr_cutoff = 0.7
alpha = 0.05

study = 'Lifeline_deep'
run_type = 'within'

base_path = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas/{study}/{run_type}'
sig_file = os.path.join(base_path, 'mb_gwas_significant.h5')
data_file = os.path.join(base_path, 'raw_data', 'mb_gwas_{X}_{Y}.h5' if run_type == 'within' else 'mb_gwas_Rep_all_{Y}.h5')
clump_file = os.path.join(base_path, 'clumping2', 'mb_gwas_{X}_{Y}.pkl')


def clump(X, Y):

    sig_df = pd.read_hdf(sig_file)
    sig_df = sig_df[(sig_df.index.get_level_values('Species') == X) & (sig_df.index.get_level_values('Y') == Y)]
    # sig_df = pd.concat([sig_df.head(100), sig_df.tail(5)])
    sig_df = sig_df.sort_values(pval_col, ascending=True)

    data = pd.read_hdf(data_file.replace('{X}', X).replace('{Y}', Y))
    data = data[(data.index.get_level_values('Species') == X)]  # makes it much more memory efficient
    data = data.reset_index(sig_df.index.names)
    data['Position'] = data['Position'].astype(int)
    data = data.groupby(sig_df.index.names)[data.columns[-1]]

    sig_snps = []
    all_snps = {}
    for new_snp in sig_df.index:
        add = True
        for sig_snp in sig_snps:
            if (new_snp[2].split('_')[1] == sig_snp[2].split('_')[1]) & (abs((new_snp[3] - sig_snp[3])) <= max_positions):
                common_samples = set(data.get_group(new_snp).index) & set(data.get_group(sig_snp).index)
                r, p = pearsonr(data.get_group(new_snp).loc[common_samples], data.get_group(sig_snp).loc[common_samples])
                if (abs(r) >= abs_corr_cutoff) & (p < alpha):
                    add = False
                    break
        if add:
            sig_snps.append(new_snp)
            all_snps[new_snp] = new_snp
        else:
            all_snps[new_snp] = sig_snp

    for i, snp2 in enumerate(sig_snps[::-1]):  # from weak to strong
        for snp1 in sig_snps[:-i-1]:  # from strong to weak
            common_samples = set(data.get_group(snp1).index) & set(data.get_group(snp2).index)
            r, p = pearsonr(data.get_group(snp1).loc[common_samples], data.get_group(snp2).loc[common_samples])
            if (abs(r) >= abs_corr_cutoff) & (p < alpha):
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

    # os.makedirs(os.path.dirname(clump_file))#, exist_ok=True)

    with qp(jobname='WLNclumping', _tryrerun=True, _mem_def='5G') as q:
        q.startpermanentrun()
        tkttores = {}

        significant_df = pd.read_hdf(sig_file)
        # significant_df = pd.concat([significant_df.head(100), significant_df.tail(5)])

        print('start sending jobs')
        runs = significant_df.reset_index()[['Species', 'Y']].drop_duplicates().reset_index(drop=True)
        for i, (speciesX, speciesY) in runs.iterrows():
            if not os.path.exists(clump_file.replace('{X}', speciesX).replace('{Y}', speciesY)):
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
        significant_df.to_hdf(sig_file.replace('.h5', '_clumping2.h5'), key='snps')
        print('finished df update')
