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

    if not kwargs['draw_pval']:
        sp_df = sp_df[sp_df[kwargs['pval_col']] < kwargs['pval_cutoff']]
    sp_df[kwargs['coef_col']] = sp_df[kwargs['coef_col']]*MAF_1_VALUE

    d = dict()

    black = (0.0, 0.0, 0.0, 0.7)
    white = (1.0, 1.0, 1.0, 0.7)
    blue = (0.0, 0.0, 1.0, 0.7)
    red = (1.0, 0.0, 0.0, 0.7)

    color_dict = {True: red, False: blue}
    color_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}

    d['marker'] = 'o'
    d['s'] = -np.log10(sp_df[kwargs['pval_col']]).values/2  # to not divide by max because then it is species subjective
    d['facecolor'] = (sp_df[kwargs['coef_col']] > 0).map(color_dict).values
    d['edgecolor'] = black
    d['linewidths'] = 1
    d['alpha'] = 0.5
    if 'text' in sp_df.columns:
        d['text'] = sp_df['text']

    if kwargs['legend_elements'] is None:
        legend_coefficient = [Line2D([0], [0], linewidth=0, label=color_label_dict[k], alpha=d['alpha'],
                                     markersize=10, marker='o', markeredgewidth=0,
                                     markerfacecolor=color_dict[k], markeredgecolor=white)
                              for k in color_dict.keys()]

        legend = legend_coefficient

    else:
        legend = kwargs['legend_elements']

    return sp_df, d, legend


def text_func_annotated_between(df: pd.DataFrame, **kwargs):
    df = df.reset_index().drop_duplicates(list(df.index.names))
    text = f"{df.shape[0]:,} SNPs in {df['Species'].unique().shape[0]:,} species\n" + \
           f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,} SNPs passed {kwargs['pval_cutoff']:.2e} p-value cutoff"
    return text


def text_func_annotated_within(df: pd.DataFrame, **kwargs):
    df = df.reset_index().drop_duplicates(list(df.index.names))
    text = f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,}/{df.shape[0]:,} SNPs passed {kwargs['pval_cutoff']:.2e} p-value cutoff"
    return text


if __name__ == '__main__':

    study = '10K'
    run_type = 'between'

    input_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas/{study}/{run_type}'
    output_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/figs/{run_type}_species'
    jobs_dir = os.path.join(input_dir, 'jobs')

    data_df = None
    annotations_df = None#pd.read_hdf(os.path.join(input_dir, 'snps_gene_annotations_short.h5'))[['text']]

    counts = pd.read_hdf(os.path.join(input_dir, 'mb_gwas_counts.h5'))
    alpha = 0.01/counts.sum()

    # queue
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    with qp(jobname=run_type, _tryrerun=True) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for file in glob.glob(os.path.join(input_dir, 'raw_hdfs', 'mb_gwas_Rep_*_Rep_*.h5')):
            kwargs = {'mwas_fname': file,
                      'out_dir': output_dir,
                      'data_fname': data_df,
                      'annotations_df': annotations_df,
                      'manhattan_draw_func': color_by_coef,
                      'manhattan_text_func': text_func_annotated_within if run_type == 'within' else
                                             text_func_annotated_between,
                      'maf_col': 'MAF',
                      'coef_col': 'Coef',
                      'pval_col': 'Pval',
                      'pval_cutoff': alpha}
            tkttores[file] = q.method(mwas_plots.run, kwargs=kwargs)
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
