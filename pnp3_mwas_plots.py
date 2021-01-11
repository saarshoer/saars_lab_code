import os
import glob
import mwas_plots
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers


def draw_func(sp_df: pd.DataFrame, **kwargs):
    d = dict()

    black = (0.0, 0.0, 0.0, 0.7)
    white = (1.0, 1.0, 1.0, 0.7)
    blue = (0.0, 0.0, 1.0, 0.7)
    red = (1.0, 0.0, 0.0, 0.7)

    facecolor_dict = {True: red, False: blue}
    facecolor_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}

    d['marker'] = 'o'
    d['s'] = -np.log10(
        sp_df[kwargs['pval_col']]).values  # to not divide by max because than it is species subjective
    # d['facecolor'] = sp_df['annotation'].map(facecolor_dict).values
    d['facecolor'] = (sp_df['Coef'] > 0).map(facecolor_dict).values
    d['edgecolor'] = black
    d['linewidths'] = 1
    d['alpha'] = 0.5
    d['text'] = sp_df['gene'].fillna('')

    if kwargs['legend_elements'] is None:
        legend_coefficient = [Line2D([0], [0], linewidth=0, label=facecolor_label_dict[k], alpha=d['alpha'],
                                     markersize=10, marker='s', markeredgewidth=0,
                                     markerfacecolor=facecolor_dict[k], markeredgecolor=white) for k in
                              facecolor_dict.keys()]

        legend = legend_coefficient
    else:
        legend = kwargs['legend_elements']

    return sp_df, d, legend


def text_func(df: pd.DataFrame, **kwargs):
    text = f"{df.shape[0]:,} SNPs in {df.index.get_level_values('Species').unique().shape[0]:,} species\n" + \
           f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,} SNPs passed {kwargs['pval_cutoff']:.2e} {kwargs['pval_col'].replace('_', ' ').replace('Pval', 'p-value')} cutoff"
    return text

run_type = 'PNP3_gut_mwas_0months_delta_species_exist'

input_path = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/all', run_type, 'mb_gwas_per_y')
output_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/figs/MWAS', run_type)
jobs_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/jobs/')
annotations_path = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/all', run_type, 'annotations_significant', 'snps_amino_acids.h5')

# annotations
print('starting annotations')
annotations_df = pd.read_hdf(annotations_path)
print('finished annotations')

# queue
os.chdir(jobs_path)
sethandlers(file_dir=jobs_path)

with fakeqp(jobname=run_type, _delete_csh_withnoerr=True, q=['himem7.q'], _mem_def='10G', _tryrerun=True, _num_reruns=3) as q:
    q.startpermanentrun()
    tkttores = {}

    print('start sending jobs')
    for y in ['bt__alt_gpt', 'bt__ast_got', 'bt__crp_hs', 'bt__crp_synthetic', 'bt__hemoglobin', 'bt__neutrophils_',
              'bt__wbc', 'carbohydrates', 'hips', 'sitting_blood_pressure_pulse_rate']:
        file = os.path.join(input_path, f'mb_gwas_{y}.h5')
        kwargs = {'mwas_fname': file, 'annotations_df': annotations_df, 'out_dir': output_path,
                  'manhattan_draw_func': draw_func,
                  'manhattan_text_func': text_func,
                  'pval_col': 'Global_Bonferroni', 'pval_cutoff': 0.05}
        tkttores[file] = q.method(mwas_plots.run, kwargs=kwargs)
    print('finished sending jobs')

    print('start waiting for jobs')
    for k, v in tkttores.items():
        q.waitforresult(v)
    print('finished waiting for jobs')

print('done')
