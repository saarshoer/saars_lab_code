import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from analysis import segal_name
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE

study = '10K'
min_periods = 50
method = 'pearson'
# min_reads_per_snp = 3

base_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/{study}'
maf_file = os.path.join(base_dir, 'pcs_covariate', 'mb_snp_maf_{}.h5')
# cov_file = os.path.join(f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/{study}/coverage_covariate', 'mb_snp_coverage_{}.h5')

res_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/strains/{study}/linkage_disequilibrium'


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
        maf_filtered = maf_filtered.where(np.triu(np.ones(maf_filtered.shape, dtype='uint8'), k=1).astype(bool))

        maf_filtered = maf_filtered.stack()
        maf_filtered.to_hdf(os.path.join(res_dir, f'{species}_C_{contig}.h5'), key='snps', complevel=9)

        # chunks = list(np.arange(0, maf_filtered.shape[1], 5000)) + [maf_filtered.shape[1]]
        # for i in np.arange(len(chunks) - 1):
        #     maf_filtered.iloc[:, chunks[i]:chunks[i+1]].stack().to_hdf(os.path.join(res_dir, f'{species}_C_{contig}.h5'), key=str(i), complevel=9)


def plot(species, files, files2, files3):

    def get(fs):

        r = []
        for file in fs:
            # contig = file.split('_')[-1].split('.')[0]
            with pd.HDFStore(file, 'r') as hdf:
                for key in hdf.keys():
                    r.append(hdf[key].to_frame('r').sample(frac=0.4, random_state=42))#.assign(species=species, contig=contig)) for memory efficiency
            del hdf
        if len(r) > 0:
            r = pd.concat(r)
            r['ld'] = r['r'] * r['r']
            r['distance'] = r.index.get_level_values(1) - r.index.get_level_values(0)

            r['ld_bin'] = pd.cut(x=r['ld'], bins=ld_bin, labels=ld_bin[:-1])
            r['distance_bin'] = pd.cut(x=np.log10(r['distance'].astype(float)).clip(upper=6),
                                             bins=distance_bin, labels=distance_bin[:-1])

            r = r[['ld', 'ld_bin', 'distance_bin']].reset_index(drop=True)  # for memory efficiency

        return r

    ld_bin = np.round(np.arange(0, 1.02, 0.01).tolist(), 2)
    distance_bin = np.round(np.arange(0, 6.2, 0.1).tolist(), 2)

    results = get(files)

    # # just for df
    # r50 = results.groupby('distance_bin')['ld'].quantile(0.5) \
    #                          .rolling(3, min_periods=1, center=True).mean() \
    #                          .to_frame(species)
    # r90 = results.groupby('distance_bin')['ld'].quantile(0.9) \
    #                          .rolling(3, min_periods=1, center=True).mean() \
    #                          .to_frame(species)
    # return r50, r90

    heatmap_data = results.groupby(['ld_bin', 'distance_bin']).apply(len)
    heatmap_data = heatmap_data.unstack('distance_bin').fillna(0).iloc[::-1]
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

    ax = sns.heatmap(heatmap_data, xticklabels=10, yticklabels=10,
                     vmax=20, cmap='Greys', cbar_kws={'label': 'Density [%]'})
    ax.set_title(f'{species}\n{segal_name(species)[0]}')
    ax.set_xlabel('Genomic distance [log10 bp]')
    ax.set_ylabel('Linkage disequilibrium [r^2]')

    del heatmap_data

    ax2 = plt.twinx().twiny()

    data = []
    for s in ['Israel', 'Netherlands', 'United States']:

        if s == 'Netherlands':
            results = get(files2)
        elif s == 'United States':
            results = get(files3)
        if len(results) == 0:
            continue

        for p in [0.5, 0.9]:
            data.append((results.groupby('distance_bin')['ld'].quantile(p) \
                        .rolling(3, min_periods=1, center=True).mean() \
                        .to_frame()*100).reset_index().assign(Color=f'LD {p*100:.0f}').assign(Line=s))
        del results
    data = pd.concat(data).reset_index(drop=True)

    sns.lineplot(x='distance_bin', y='ld', hue='Color', style='Line', data=data, ax=ax2,
                 hue_order=[f'LD 50', 'LD 90'], style_order=['Israel', 'Netherlands', 'United States'])

    ax2.legend(loc='upper right')
    ax2.set_xlabel('')
    ax2.set_xticks([])
    ax2.set_xlim(0.2, 6)
    ax2.set_ylabel('')
    ax2.set_yticks([])
    ax2.set_ylim(0, 100)

    plt.savefig(os.path.join(os.path.dirname(res_dir), 'figs', 'linkage_disequilibrium', species))
    plt.close()


# if __name__ == '__main__':
#
#     jobs_dir = os.path.join(os.path.dirname(res_dir), 'jobs')
#     os.chdir(jobs_dir)
#     sethandlers(file_dir=jobs_dir)
#
#     os.makedirs(res_dir, exist_ok=True)
#
#     with qp(jobname='LK', _tryrerun=True, _delete_csh_withnoerr=True, _mem_def='150G', max_r=125) as q:
#         q.startpermanentrun()
#         tkttores = {}
#
#         sc = pd.read_hdf(os.path.join(base_dir, 'within', 'mb_gwas_counts.h5'))
#         sc['Contig'] = sc.index.get_level_values('Contig').str.split('_').str[1]
#         sc = sc['Contig'].reset_index('Species').drop_duplicates().values
#
#         print('start sending jobs')
#         for s, c in sc:
#             if not os.path.exists(os.path.join(res_dir, f'{s}_C_{c}.h5')):
#                 print(s, c)
#                 tkttores[f'{s}_{c}'] = q.method(do, [s, c], _job_name=f'lk{s.split("_")[-1]}')
#         print('finished sending jobs')
#
#         print('start waiting for jobs')
#         for k, v in tkttores.items():
#             q.waitforresult(v)
#         print('finished waiting for jobs')

if __name__ == '__main__':

    jobs_dir = os.path.join(os.path.dirname(res_dir), 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(res_dir, exist_ok=True)

    with qp(jobname='LK', _tryrerun=True, _num_reruns=10, _delete_csh_withnoerr=True, max_r=125) as q:#_mem_def='100G',
        q.startpermanentrun()
        tkttores = {}

        k10_sc = pd.read_hdf(os.path.join(base_dir, 'within', 'mb_gwas_counts.h5'))
        k10_sc['Contig'] = k10_sc.index.get_level_values('Contig').str.split('_').str[1]
        k10_sc = k10_sc['Contig'].reset_index('Species').drop_duplicates()

        lld_sc = pd.read_hdf(os.path.join(base_dir.replace('10K', 'Lifeline_deep'), 'within', 'mb_gwas_counts.h5'))
        lld_sc['Contig'] = lld_sc.index.get_level_values('Contig').str.split('_').str[1]
        lld_sc = lld_sc['Contig'].reset_index('Species').drop_duplicates()

        d2u_sc = pd.read_hdf(os.path.join(base_dir.replace('10K', 'D2'), 'within', 'mb_gwas_counts.h5'))
        d2u_sc['Contig'] = d2u_sc.index.get_level_values('Contig').str.split('_').str[1]
        d2u_sc = d2u_sc['Contig'].reset_index('Species').drop_duplicates()

        # results50 = None#pd.read_pickle(os.path.join(os.path.dirname(res_dir), 'data_frames', 'linkage_disequilibrium_50.df'))
        # results90 = None#pd.read_pickle(os.path.join(os.path.dirname(res_dir), 'data_frames', 'linkage_disequilibrium_90.df'))

        print('start sending jobs')
        for s, cs in k10_sc.groupby('Species')['Contig']:
            # if s in results90.columns:
            #     continue
            # if not os.path.exists(os.path.join(os.path.dirname(res_dir), 'figs', 'linkage_disequilibrium', f'{s}.png')):
            k10_f = [os.path.join(res_dir, f'{s}_C_{c}.h5') for c in cs]
            lld_f = [os.path.join(res_dir.replace('10K', 'Lifeline_deep'), f'{s}_C_{c}.h5') for c in
                     lld_sc.loc[lld_sc['Species'] == s, 'Contig']]
            d2u_f = [os.path.join(res_dir.replace('10K', 'D2'), f'{s}_C_{c}.h5') for c in
                     d2u_sc.loc[d2u_sc['Species'] == s, 'Contig']]
            # print(s)
            tkttores[s] = q.method(plot, [s, k10_f, lld_f, d2u_f], _job_name=f'lk{s.split("_")[-1]}')
        del k10_sc, lld_sc, d2u_sc
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            try:
                q.waitforresult(v, _assert_on_errors=False)
                # results = q.waitforresult(v, _assert_on_errors=False)
                # results50 = results50.join(results[0], how='outer') if results50 is not None else results[0]
                # results90 = results90.join(results[1], how='outer') if results90 is not None else results[1]
            except:
                print(k)
                continue
        # results50.to_pickle(os.path.join(os.path.dirname(res_dir), 'data_frames', 'linkage_disequilibrium_50.df'))
        # results90.to_pickle(os.path.join(os.path.dirname(res_dir), 'data_frames', 'linkage_disequilibrium_90.df'))
        print('finished waiting for jobs')
