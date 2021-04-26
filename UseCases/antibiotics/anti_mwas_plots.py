import os
import glob
import mwas_plots
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers


def color_by_coef(sp_df: pd.DataFrame, **kwargs):
    d = dict()

    black = (0.0, 0.0, 0.0, 0.7)
    white = (1.0, 1.0, 1.0, 0.7)
    blue = (0.0, 0.0, 1.0, 0.7)
    red = (1.0, 0.0, 0.0, 0.7)

    color_dict = {True: red, False: blue}
    color_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}

    d['marker'] = 'o'
    d['s'] = -np.log10(sp_df[kwargs['pval_col']]).values/2  # to not divide by max because than it is species subjective
    d['facecolor'] = (sp_df['Coef'] > 0).map(color_dict).values
    d['edgecolor'] = black
    d['linewidths'] = 1
    d['alpha'] = 0.5

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
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,} SNPs passed {kwargs['pval_cutoff']:.2e} {kwargs['pval_col'].replace('_', ' ').replace('Pval', 'p-value')} cutoff"
    return text


def text_func_annotated_within(df: pd.DataFrame, **kwargs):
    df = df.reset_index().drop_duplicates(list(df.index.names))
    text = f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,}/{df.shape[0]:,} SNPs passed {kwargs['pval_cutoff']:.2e} {kwargs['pval_col'].replace('_', ' ').replace('Pval', 'p-value')} cutoff"
    return text


run_type = 'between_species'

input_dir = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed', run_type)
output_dir = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs/new', run_type)
jobs_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs/'

# queue
os.chdir(jobs_dir)
sethandlers(file_dir=jobs_dir)

with qp(jobname=run_type, _tryrerun=True) as q:
    q.startpermanentrun()
    tkttores = {}

    print('start sending jobs')
    for file in glob.glob(os.path.join(input_dir, 'SGB_*.h5')):  # 9710, 10068
        kwargs = {'mwas_fname': file, 'out_dir': output_dir,
                  'manhattan_draw_func': color_by_coef,
                  'manhattan_text_func': text_func_annotated_between if run_type == 'between_species' else text_func_annotated_within,
                  'pval_col': 'Pval', 'pval_cutoff': 0.05/26068850133}
        tkttores[file] = q.method(mwas_plots.run, kwargs=kwargs)
    print('finished sending jobs')

    print('start waiting for jobs')
    for k, v in tkttores.items():
        q.waitforresult(v)
    print('finished waiting for jobs')

print('done')
