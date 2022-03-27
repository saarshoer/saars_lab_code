import os
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from sklearn.decomposition import PCA
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MBSNPLoader, get_mbsnp_loader_class, MAF_1_VALUE, MAF_MISSING_VALUE

df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/data_frames'
pca_output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_pca'
cor_output_dir = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_corr'


def _filter_snps_by_stats_of_major(x):
    _min_on_minor_per_snp = 50
    _min_subjects_per_snp = 500

    major_count = MBSNPLoader.compute_maf_major_counts(x).values
    maf_counts = MBSNPLoader.compute_maf_counts(x).values
    res = pd.Series(True, index=x.columns)
    if _min_on_minor_per_snp is not None:
        res.iloc[:] = np.logical_and(res.iloc[:],
                                     np.logical_and(maf_counts - major_count >= _min_on_minor_per_snp,
                                                    major_count >= _min_on_minor_per_snp))
    if _min_subjects_per_snp is not None:
        res.iloc[:] = np.logical_and(res.iloc[:], maf_counts >= _min_subjects_per_snp)
    x = x[res[res].index]
    return x


def do(s):

    samples_set = pd.read_pickle(os.path.join(df_dir, 'snps.df')).index.tolist()
    snp_loader = get_mbsnp_loader_class('Gut')

    all_contigs = None
    for contig_dl in \
            snp_loader.get_data(df=s,
                                min_reads_per_snp=1,
                                min_positions_per_sample=1,
                                min_samples_per_snp=500,
                                min_samples_per_snp_cached=500,
                                data_gen_f='_species_maf_contig_iter',
                                max_on_fraq_major_per_snp=0.8,
                                largest_sample_per_user=False,
                                subjects_df=None,
                                groupby_reg='first',
                                snp_set=None,
                                samples_set=samples_set,
                                filter_by_species_existence=False,
                                subsample_dir='10K',
                                column_clusterer=None):

        contig_dl.df.columns = [f'{contig_dl.added_data}_pos{c}' for c in contig_dl.df.columns]
        if all_contigs is None:
            all_contigs = contig_dl.df
        else:
            all_contigs = all_contigs.join(contig_dl.df)

    all_contigs = all_contigs.fillna(MAF_MISSING_VALUE)
    all_contigs = _filter_snps_by_stats_of_major(all_contigs)

    # pca
    pca_all_contigs = all_contigs.replace(MAF_MISSING_VALUE, MAF_1_VALUE)####maybe not replace so missing could be reflected by pca, or use more complex imputation
    pca = PCA(n_components=10)
    pd.DataFrame(pca.fit_transform(pca_all_contigs),
                 index=pca_all_contigs.index, columns=[f'PC{i}' for i in np.arange(1, 11)]) \
        .to_pickle(os.path.join(pca_output_dir, f'{s}.df'))
    pd.DataFrame(pca.explained_variance_)\
        .to_csv(os.path.join(pca_output_dir, f'{s}.csv'))
    del pca_all_contigs

    # correlation
    cor_all_contigs = all_contigs.replace(MAF_MISSING_VALUE, np.nan)
    cor_all_contigs = cor_all_contigs.corr()###min_periods=500
    cor_all_contigs = cor_all_contigs.where(pd.np.triu(pd.np.ones(cor_all_contigs.shape[1]), k=1).astype(bool)).stack()
    cor_all_contigs = cor_all_contigs[cor_all_contigs.abs() > 0.8]
    cor_all_contigs.to_hdf(os.path.join(cor_output_dir, f'{s}.h5'), key='pos')
    del cor_all_contigs


if __name__ == '__main__':
    jobs_dir = os.path.join(pca_output_dir, 'jobs')
    os.makedirs(jobs_dir, exist_ok=True)
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    species_set = ['Rep_959']#pd.read_pickle(os.path.join(df_dir, 'snps.df')).columns.tolist()

    with qp(jobname='ps', max_r=20, _mem_def='20G', _tryrerun=True) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for species in species_set:
            tkttores[species] = q.method(do, [species], _job_name=f'ps{species.split("_")[-1]}')
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')
