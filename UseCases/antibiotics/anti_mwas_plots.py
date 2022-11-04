import os
import glob
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_plots
from LabData.DataLoaders.MBSNPLoader import MAF_1_VALUE


def color_by_coef(sp_df: pd.DataFrame, **kwargs):

    if kwargs['draw_col'] != kwargs['pval_col']:
        sp_df = sp_df[sp_df[kwargs['pval_col']] < kwargs['pval_cutoff']]
        if 'coef' in kwargs['draw_col'].lower():
            sp_df[kwargs['coef_col']] = sp_df[kwargs['coef_col']]*MAF_1_VALUE

    d = dict()

    blue = (0.0, 0.0, 1.0)
    red = (1.0, 0.0, 0.0)

    color_dict = {True: red, False: blue}
    color_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}

    d['marker'] = 'o'
    d['s'] = -np.log10(sp_df[kwargs['pval_col']]).values/2  # to not divide by max because then it is species subjective

    d['edgecolor'] = 'black'
    if 'validation_level' in sp_df.columns:
        d['linewidths'] = (sp_df['validation_level'] == 'CI overlap').replace({True: 1, False: 0})
    else:
        d['linewidths'] = 1

    if 'clumping' in sp_df.columns:
        is_clumping_rep = sp_df.index == sp_df['clumping']
        d['facecolor'] = (sp_df['Coef'] > 0).map(color_dict)
        d['facecolor'].loc[is_clumping_rep] = [(c[0], c[1], c[2], 0.8) for c in d['facecolor'].loc[is_clumping_rep]]
        d['facecolor'].loc[~is_clumping_rep] = [(c[0], c[1], c[2], 0.2) for c in d['facecolor'].loc[~is_clumping_rep]]
        d['facecolor'] = d['facecolor'].values
        d['alpha'] = None
        if 'text' in sp_df.columns:
            sp_df.loc[~is_clumping_rep, 'text'] = None
    else:
        d['facecolor'] = (sp_df['Coef'] > 0).map(color_dict).values
        d['alpha'] = 0.5
    if 'text' in sp_df.columns:
        d['text'] = sp_df['text']

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

    text = f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,}/{df.shape[0]:,} " + \
           f"SNPs passed {kwargs['pval_cutoff']:.2e} p-value cutoff"
    if 'validation_level' in df.columns:
        text = text + f"\n{df[df['validation_level'] == 'CI overlap'].shape[0]:,}/" + \
                        f"{df[df['validation_level'] != 'not tested'].shape[0]:,} " + \
                        f"SNPs CI overlaps with validation cohort"
    if 'clumping' in df.columns:
        text = text + f"\n{df['clumping'].unique().shape[0]:,} clumped SNPs"

    if (df['Y'] != df['Species']).any():
        text = f"{df.shape[0]:,} SNPs in {df['Species'].unique().shape[0]:,} species\n" + text

    return text


if __name__ == '__main__':

    study = '10K'
    run_type = 'within'

    input_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas/{study}/{run_type}'
    output_dir = os.path.join(input_dir, 'figs', 'species')
    jobs_dir = os.path.join(input_dir, 'jobs')

    ydata_fname = os.path.join(os.path.dirname(input_dir), 'data_frames', 'abundance.df')

    mwas_df = pd.read_hdf(os.path.join(input_dir, 'mb_gwas_significant_full.h5'))[['clumping', 'validation_level']]
    annotations_df = pd.read_hdf(os.path.join(input_dir, 'snps_gene_annotations_short.h5'))[['text']]
    annotations_df = mwas_df.join(annotations_df)
    del mwas_df

    counts = pd.read_hdf(os.path.join(input_dir, 'mb_gwas_counts.h5'))
    alpha = 0.01/counts.sum().iloc[0]
    del counts

    fmt = lambda coef, pos: f'$2^{{{coef:.1f}}}$'

    # queue
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    with qp(jobname=run_type, _tryrerun=True) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for mwas_fname in glob.glob(os.path.join(input_dir, 'raw_hdfs', 'mb_gwas_Rep_*_Rep_*.h5')):
            data_fname = mwas_fname.replace('raw_hdfs', 'raw_data')
            kwargs = {'mwas_fname': mwas_fname,
                      'data_fname': data_fname if os.path.exists(data_fname) else None,
                      'ydata_fname': ydata_fname,
                      'out_dir': output_dir,
                      'annotations_df': annotations_df,
                      'manhattan_draw_func': color_by_coef,
                      'manhattan_text_func': text_func,
                      'r2_col': 'R2',
                      'maf_col': 'MAF',
                      'coef_col': 'Coef',
                      'pval_col': 'Pval',
                      'pval_cutoff': alpha,
                      'fmt': fmt}
            tkttores[mwas_fname] = q.method(mwas_plots.run, kwargs=kwargs)
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
