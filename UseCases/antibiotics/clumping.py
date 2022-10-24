import os
import pickle
import pandas as pd
from scipy.stats import pearsonr
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

pval_col = 'Pval'
# method = 'pearson'
corr_cutoff = 0.75
min_periods = 250

study = '10K'
run_type = 'within'

base_path = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas/{study}/{run_type}'
sig_file = os.path.join(base_path, 'mb_gwas_significant.h5')
data_file = os.path.join(base_path, 'raw_data', 'mb_gwas_{X}_{Y}.h5')
clump_file = os.path.join(base_path, 'clumping', 'mb_gwas_{X}_{Y}.pkl')


def clump(X, Y):

    def long_name(snp):
        return Y, X, snp[0], snp[1]

    sig_df = pd.read_hdf(sig_file)
    # sig_df = pd.concat([sig_df.head(5), sig_df.tail(5)])
    sig_df = sig_df.xs(Y, level='Y').xs(X, level='Species')
    sig_df = sig_df.sort_values(pval_col, ascending=True)

    data_df = pd.read_hdf(data_file.replace('{X}', X).replace('{Y}', Y))
    data_df = data_df.xs(Y, level='Y').xs(X, level='Species')

    sig_snp = None
    sig_snps = []
    all_snps = {}
    for new_snp in sig_df.index:
        new_snp_data = data_df.xs(new_snp[0], level='Contig').xs(str(new_snp[1]), level='Position')
        add = True
        for sig_snp in sig_snps:
            sig_snp_data = data_df.xs(sig_snp[0], level='Contig').xs(str(sig_snp[1]), level='Position')
            # corr = new_snp_data.corr(sig_snp_data, method=method, min_periods=min_periods)  # too slow
            shared_samples = new_snp_data.index.intersection(sig_snp_data.index)
            if len(shared_samples) >= min_periods:
                corr, pval = pearsonr(new_snp_data.loc[shared_samples], sig_snp_data.loc[shared_samples])
            else:
                corr, pval = 0, 1
            if corr >= corr_cutoff:
                add = False
                break
        if add:
            sig_snps.append(new_snp)
        all_snps[long_name(new_snp)] = long_name(new_snp) if add else long_name(sig_snp)

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

    with qp(jobname='clumping', _tryrerun=True, _specific_nodes='plink') as q:
        q.startpermanentrun()
        tkttores = {}

        significant_df = pd.read_hdf(sig_file)
        # significant_df = pd.concat([significant_df.head(5), significant_df.tail(5)])

        print('start sending jobs')
        runs = significant_df.reset_index()[['Species', 'Y']].drop_duplicates().reset_index(drop=True)
        for i, (speciesX, speciesY) in runs.iterrows():
            tkttores[i] = q.method(clump, (speciesX, speciesY))
        print('finished sending jobs')

        print('start waiting for jobs')
        results_files = []
        for k, v in tkttores.items():
            results_files.append(q.waitforresult(v))
        print('finished waiting for jobs')
        # results_files = glob.glob(clump_file.replace('{X}', '*').replace('{Y}', '*'))

        print('start df update')
        significant_df['clumping'] = None
        for results_file in results_files:
            with open(results_file, 'rb') as f:
                run_results = pickle.load(f)
            significant_df.loc[run_results.keys()] = significant_df.loc[run_results.keys()].assign(
                clumping=run_results.values())
        significant_df.to_hdf(sig_file.replace('.h5', '_clumping.h5'), key='snps')
        print('finished df update')
