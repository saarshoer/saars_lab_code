import os
import glob
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_plots
from LabData.DataLoaders.MBSNPLoader import MAF_1_VALUE

all_contigs = pd.read_csv('/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/Unicorn/URS/URS_Build/Rep_all/all_contigs.csv')
all_contigs['Species'] = 'Rep_' + all_contigs['contig'].str.split('_').str[1]
all_contigs['contig_id'] = all_contigs['contig'].str.split('_').str[-1].astype(int)
all_contigs = all_contigs.set_index(['Species', 'contig_id'])['len'].sort_values().to_frame('contig_len')*-1

# dna_genes = pd.read_csv('/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/dna_genes.csv')
# dna_genes['Species'] = 'Rep_' + dna_genes['Species'].astype(str)
# dna_genes = dna_genes.set_index(['Species', 'contig_id', 'start_pos'])
#
# mobilome = pd.read_csv('/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/mobilome.csv')
# mobilome = mobilome.set_index(['GeneID'])


def color_by_coef(sp_df: pd.DataFrame, **kwargs):
    # sp_df = sp_df.join(dna_genes[dna_genes.index.get_level_values('Species').isin(sp_df.index.get_level_values('Species'))],
    #                    on=['Species', 'contig_id', 'Position'], how='outer')  # before contig to be ordered
    sp_df = sp_df.join(all_contigs, on=['Species', 'contig_id']).sort_values(['contig_len', 'contig_id', 'Position'])
    sp_df = sp_df.reset_index().drop_duplicates().set_index(sp_df.index.names)

    if kwargs['draw_col'] != kwargs['pval_col']:
        sp_df = sp_df[(sp_df[kwargs['pval_col']] < kwargs['pval_cutoff'])]# | (~sp_df['dna'].isna())
        if 'coef' in kwargs['draw_col'].lower():
            sp_df[kwargs['coef_col']] = sp_df[kwargs['coef_col']]*MAF_1_VALUE

    d = dict()

    blue = (0.0, 0.0, 1.0)
    red = (1.0, 0.0, 0.0)

    color_dict = {True: blue, False: red}
    color_label_dict = {True: 'Positive effect', False: 'Negative effect'}

    d['marker'] = 'o'
    d['s'] = -np.log10(sp_df[kwargs['pval_col']]).values/2  # to not divide by max because then it is species subjective

    d['edgecolor'] = 'black'
    if 'validation_level' in sp_df.columns:
        d['linewidths'] = (sp_df['validation_level'] == 'CI overlap').replace({True: 1, False: 0})
    else:
        d['linewidths'] = 0

    d['facecolor'] = (sp_df['Coef'] > 0).map(color_dict)
    if 'clumping' in sp_df.columns:
        is_clumping_rep = sp_df.index == sp_df['clumping']
        d['facecolor'].loc[is_clumping_rep] = [(c[0], c[1], c[2], 0.8) for c in d['facecolor'].loc[is_clumping_rep]]
        d['facecolor'].loc[~is_clumping_rep] = [(c[0], c[1], c[2], 0.2) for c in d['facecolor'].loc[~is_clumping_rep]]
        d['alpha'] = None
        # if 'text' in sp_df.columns:
        #     sp_df.loc[~is_clumping_rep, 'text'] = None
    elif 'feature_importance' in sp_df.columns:
        is_important = sp_df['feature_importance'] > 0
        d['facecolor'].loc[is_important] = [(c[0], c[1], c[2], 0.8) for c in d['facecolor'].loc[is_important]]
        d['facecolor'].loc[~is_important] = [(c[0], c[1], c[2], 0.2) for c in d['facecolor'].loc[~is_important]]
        d['alpha'] = None
    else:
        d['alpha'] = 0.5
    d['facecolor'] = d['facecolor'].values

    if 'text' in sp_df.columns:
        d['text'] = sp_df['text']

    # if (sp_df.index.get_level_values('Y') == sp_df.index.get_level_values('Species')).any():
    #
    #     sp_df['tick'] = '|contig start'
    #     # sp_df['tick'] = '|contig start\ndna gene\nmobilome' if 'GeneID' in sp_df.columns else '|contig start\ndna gene'
    #
    #     # contigs
    #     sp_df['tick'][1:] = sp_df['contig_id'][1:].values != sp_df['contig_id'][:-1].values
    #     sp_df['tick'] = sp_df['tick'].replace({True: '|', False: ''})#{True: '|\n', False: '\n'}
    #
    #     # # dna genes
    #     # sp_df['tick'] = sp_df['tick'] + sp_df['dna'].fillna('') + '\n'
    #
    #     # # mobilome
    #     # if 'GeneID' in sp_df.columns:
    #     #     sp_df = sp_df.join(mobilome, on='GeneID')
    #     #     sp_df['tick'] = sp_df['tick'] + sp_df['mobilome'].fillna('')
    #
    #     # sp_df['tick'] = sp_df['tick'].replace({'\n\n': None})
    #     d['tick'] = sp_df['tick']

    if kwargs['legend_elements'] is None:
        legend_coefficient = [Line2D([0], [0], linewidth=0, label=color_label_dict[k], alpha=d['alpha'],
                                     markersize=10, marker='o', markeredgewidth=0,
                                     markerfacecolor=color_dict[k], markeredgecolor='white')
                              for k in color_dict.keys()]

        legend = legend_coefficient

    else:
        legend = kwargs['legend_elements']

    return sp_df, d, legend


def text_func(df: pd.DataFrame, **kwargs):
    df = df.reset_index().drop_duplicates(list(df.index.names))
    is_sig = df[kwargs['pval_col']] <= kwargs['pval_cutoff']
    n_sig = f"{sum(is_sig):,}"

    text = ""

    if (df['Y'] != df['Species']).any():
        text = f"{df['Species'].unique().shape[0]:,} species\n" + text

    text = text + f"{df['N'].mean():,.0f}±{df['N'].std():,.0f} samples per SNP\n"\

    text = text + f"{n_sig}/{df.shape[0]:,} significant SNPs\n"

    if 'clumping' in df.columns and n_sig != '0':
        text = text + f"{df['clumping'].unique().shape[0]:,}/{n_sig} independent SNPs\n"

    if 'validation_level' in df.columns:
        text = text + f"{df[is_sig & (df['validation_level'] == 'CI overlap')].shape[0]:,}/" + \
                      f"{df[is_sig & (df['validation_level'] != 'not tested')].shape[0]:,} " + \
                        f"validated SNPs"

    return text


if __name__ == '__main__':

    study = '10K'
    run_type = 'within'
    data_plots = True

    input_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/{study}/{run_type}'
    output_dir = os.path.join(input_dir, 'figs_IHMC', 'species')
    jobs_dir = os.path.join(input_dir, 'jobs')

    ydata_fname = os.path.join(os.path.dirname(input_dir), 'data_frames', 'abundance.df')

    mwas_df = pd.read_hdf(os.path.join(input_dir, 'mb_gwas_significant_clumping_validation.h5')).assign(validation_level=None)
    mwas_df.loc[(mwas_df[f'validation_level_Lifeline_deep'] == 'CI overlap') & (mwas_df[f'validation_level_D2'] == 'CI overlap'), 'validation_level'] = 'CI overlap'
    mwas_df.loc[(mwas_df[f'validation_level_Lifeline_deep'] == 'not tested') | (mwas_df[f'validation_level_D2'] == 'not tested'), 'validation_level'] = 'not tested'
    mwas_df = mwas_df[['Pval', 'clumping', 'validation_level']]
    mwas_df = mwas_df[mwas_df.index.get_level_values('Species') == 'Rep_620']
    annotations_df = pd.read_hdf(os.path.join(input_dir, 'snps_gene_annotations_short.h5'))[['GeneID', 'text']]
    mwas_df = mwas_df.join(annotations_df)
    del annotations_df

    if True:#not data_plots:
        group_text = pd.DataFrame(mwas_df.sort_values('Pval').groupby(['Y', 'Species', 'GeneID']).apply(
            lambda g: g.reset_index()[list(mwas_df.index.names) + ['GeneID']].iloc[0].tolist() + [
                f'{g["text"].fillna("").iloc[0]}'+(f'({g.shape[0]})' if g.shape[0] > 1 else '')]).tolist(),
                    columns=list(mwas_df.index.names) + ['GeneID', 'text']).set_index(mwas_df.index.names)
        group_text = group_text[group_text['text'].str.contains('xerC|lptG')]
        group_text = group_text[~group_text['text'].str.split("(").str[0].duplicated()]
        group_text = group_text['text'].str.split("(").str[0]
        mwas_df = mwas_df.drop(['GeneID', 'text'], axis=1).join(group_text)
    mwas_df = mwas_df.drop(['Pval'], axis=1)

    counts = pd.read_hdf(os.path.join(input_dir, 'mb_gwas_counts.h5'))
    alpha = 0.01/counts.sum().iloc[0]
    del counts

    fmt0d = lambda coef, pos: f'$2^{{{coef:.0f}}}$'
    fmt1d = lambda coef, pos: f'$2^{{{coef:.1f}}}$'

    # queue
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    with qp(jobname=run_type, _tryrerun=True) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for mwas_fname in glob.glob(os.path.join(input_dir, 'raw_hdfs', 'mb_gwas_Rep_620_Rep_620.h5')):
            data_fname = mwas_fname.replace('raw_hdfs', 'raw_data')
            if data_plots & ~os.path.exists(data_fname):
                continue
            kwargs = {'mwas_fname': mwas_fname if not data_plots else None,
                      'data_fname': data_fname if data_plots else None,
                      'ydata_fname': ydata_fname,
                      'out_dir': output_dir,
                      'annotations_df': mwas_df,
                      'manhattan_draw_func': color_by_coef,
                      'manhattan_text_func': text_func,
                      'r2_col': 'R2',
                      'maf_col': 'MAF',
                      'coef_col': 'Coef',
                      'pval_col': 'Pval',
                      'pval_cutoff': alpha,
                      'fmt': fmt0d}
            tkttores[mwas_fname] = q.method(mwas_plots.run, kwargs=kwargs)
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
