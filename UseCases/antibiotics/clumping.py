import os
import glob
import pickle
import pandas as pd
from scipy.stats import pearsonr
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

pval_col = 'Pval'
method = 'pearson'
min_samples = 250
abs_corr_cutoff = 0.75

study = '10K'
run_type = 'between'

base_path = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas/{study}/{run_type}'
sig_file = os.path.join(base_path, 'mb_gwas_significant.h5')
data_file = os.path.join(base_path, 'raw_data', 'mb_gwas_Rep_all_{Y}.h5')
clump_file = os.path.join(base_path, 'abs_clumping', 'mb_gwas_{X}_{Y}.pkl')


def clump(X, Y):

    sig_df = pd.read_hdf(sig_file)
    sig_df = sig_df[(sig_df.index.get_level_values('Species') == X) & (sig_df.index.get_level_values('Y') == Y)]
    # sig_df = pd.concat([sig_df.head(50), sig_df.tail(5)])
    sig_df = sig_df.sort_values(pval_col, ascending=True)

    data = pd.read_hdf(data_file.replace('{X}', X).replace('{Y}', Y))
    data = data[(data.index.get_level_values('Species') == X)]  # makes it much more memory efficient
    data = data.reset_index(sig_df.index.names)
    data['Position'] = data['Position'].astype(int)
    data = data.groupby(sig_df.index.names)[[data.columns[-1]]]

    sig_snps = []
    all_snps = {}
    for new_snp in sig_df.index:
        add = True
        for sig_snp in sig_snps:
            common_samples = set(data.get_group(new_snp).index) & set(data.get_group(sig_snp).index)
            if len(common_samples) >= min_samples:
                corr = data.get_group(new_snp).corrwith(data.get_group(sig_snp), method=method)[0]
                if abs(corr) >= abs_corr_cutoff:
                    add = False
                    break
        if add:
            sig_snps.append(new_snp)
            all_snps[new_snp] = new_snp
        else:
            all_snps[new_snp] = sig_snp

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

    with qp(jobname='BKclumping', _tryrerun=True, _mem_def='5G') as q:
        q.startpermanentrun()
        tkttores = {}

        significant_df = pd.read_hdf(sig_file)
        # significant_df = pd.concat([significant_df.head(50), significant_df.tail(5)])

        print('start sending jobs')
        runs = significant_df.reset_index()[['Species', 'Y']].drop_duplicates().reset_index(drop=True)
        for i, (speciesX, speciesY) in runs.iterrows():
            # if not os.path.exists(clump_file.replace('{X}', speciesX).replace('{Y}', speciesY)):
            tkttores[i] = q.method(clump, (speciesX, speciesY))
        print('finished sending jobs')

        print('start waiting for jobs')
        results_files = []
        for k, v in tkttores.items():
            results_files.append(q.waitforresult(v))
        print('finished waiting for jobs')
        results_files = glob.glob(clump_file.replace('{X}', '*').replace('{Y}', '*'))

        print('start df update')
        significant_df['clumping'] = None
        for results_file in results_files:
            with open(results_file, 'rb') as f:
                run_results = pickle.load(f)
            significant_df.loc[run_results.keys()] = significant_df.loc[run_results.keys()].assign(
                clumping=run_results.values())
        significant_df.to_hdf(sig_file.replace('.h5', '_abs_clumping.h5'), key='snps')
        print('finished df update')
