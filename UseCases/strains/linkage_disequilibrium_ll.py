import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE

study = 'Lifeline_deep'
min_periods = 50
method = 'pearson'
# min_reads_per_snp = 3

base_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}'
maf_file = os.path.join(base_dir, 'pcs_covariate', 'mb_snp_maf_{}.h5')
# cov_file = os.path.join(f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/coverage_covariate', 'mb_snp_coverage_{}.h5')

res_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/strains/{study}/linkage_disequilibrium'


def do(species, contig=None):
    maf = pd.read_hdf(maf_file.format(species))  # contains MAF_MISSING_VALUE
    contigs = set(maf.columns.str.split('_').str[2]) if contig is None else [contig]
    for contig in contigs:
        contig_columns = maf.columns[maf.columns.str.split('_').str[2] == contig]
        contig_positions = contig_columns.str.split('_').str[0].astype(int)

        # cov = pd.read_hdf(cov_file.format(species), key=f'/C_{contig}')[contig_positions]

        # maf_filtered = maf.loc[cov.index, contig_columns]
        maf_filtered = maf[contig_columns]
        if contig == contigs[-1]:
            del maf
        maf_filtered = maf_filtered.replace(MAF_MISSING_VALUE, np.nan)
        maf_filtered.columns = contig_positions  # necessary for mask function
        # maf_filtered = maf_filtered.where(cov < min_reads_per_snp)
        maf_filtered = maf_filtered[sorted(contig_positions)]  # so pos1 will be smaller than pos2
        maf_filtered = maf_filtered.dropna(how='all', axis=0).dropna(how='all', axis=1)

        maf_filtered = maf_filtered.corr(method=method, min_periods=min_periods)
        maf_filtered = maf_filtered.where(np.triu(np.ones(maf_filtered.shape), k=1).astype(bool)).stack()
        maf_filtered.to_hdf(os.path.join(res_dir, f'{species}_C_{contig}.h5'), key='snps', complevel=9)


def plot(species, files):

    results = []
    for file in files:
        contig = file.split('_')[-1].split('.')[0]
        with open(file, 'rb') as f:
            r = pickle.load(f)
        results.append(
            pd.DataFrame(r, columns=['position', 'distance', 'n_samples', 'r', 'p']).assign(species=species, contig=contig))

    results = pd.concat(results)
    results = results[results['p'] < 0.05]
    results['ld'] = results['r'] * results['r']

    ld_bin = np.round(np.arange(0, 1.01, 0.01).tolist(), 2)
    distance_bin = np.round(np.arange(0, 5.1, 0.1).tolist(), 2)
    results['ld_bin'] = pd.cut(x=results['ld'], bins=ld_bin, labels=ld_bin[:-1])
    results['distance_bin'] = pd.cut(x=np.log10(results['distance'].astype(float)).clip(upper=5),
                                     bins=distance_bin, labels=distance_bin[:-1])

    heatmap_data = results.groupby(['ld_bin', 'distance_bin']).apply(len).unstack('distance_bin').fillna(0).iloc[::-1]
    heatmap_data = heatmap_data / heatmap_data.sum() * 100

    missing_index = list(set(ld_bin[:-1]) - set(heatmap_data.index))
    if len(missing_index) > 0:
        heatmap_data = heatmap_data.T
        heatmap_data.loc[:, missing_index] = 0
        heatmap_data = heatmap_data.T

    missing_columns = list(set(distance_bin[:-1]) - set(heatmap_data.columns))
    if len(missing_columns) > 0:
        heatmap_data.loc[:, missing_columns] = 0

    heatmap_data = heatmap_data.sort_index().T.sort_index().T.iloc[::-1]

    ax = sns.heatmap(heatmap_data, vmax=25, cmap='Greys')
    ax.set_title(species)
    ax.set_xlabel('Genomic distance [log10 bp]')
    ax.set_ylabel('LD [r^2]')

    ax2 = plt.twinx()

    for q, p in [(0.9, 'Blues'), (0.5, 'Reds')]:
        sns.lineplot(data=results.groupby('distance_bin')['ld'].quantile(q).to_frame(f'LD {q}').reset_index(),
                     palette=p, ax=ax2)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])

    # plt.show()
    plt.savefig(f'/home/saarsh/Analysis/strains/10K/figs/linkage_disequilibrium/{species}')
    plt.close()


if __name__ == '__main__':

    jobs_dir = os.path.join(os.path.dirname(res_dir), 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(res_dir, exist_ok=True)

    with qp(jobname='LL', _tryrerun=True, _delete_csh_withnoerr=True, _mem_def='10G', max_r=100) as q:
        q.startpermanentrun()
        tkttores = {}

        sc = pd.read_hdf(os.path.join(base_dir, 'within', 'mb_gwas_counts.h5'))
        sc['Contig'] = sc.index.get_level_values('Contig').str.split('_').str[1]
        sc = sc['Contig'].reset_index('Species').drop_duplicates().values

        print('start sending jobs')
        for s, c in sc:
            tkttores[f'{s}_{c}'] = q.method(do, [s, c], _job_name=f'll{s.split("_")[-1]}')
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')
