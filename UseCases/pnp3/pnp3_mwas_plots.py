import os
import mwas_plots
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


def draw_func(sp_df: pd.DataFrame, **kwargs):
    d = dict()

    black = (0.0, 0.0, 0.0, 0.7)
    white = (1.0, 1.0, 1.0, 0.7)
    blue = (0.0, 0.0, 1.0, 0.7)
    red = (1.0, 0.0, 0.0, 0.7)

    marker_dict = {'synonymous': 'o', 'non-synonymous': 'o', 'non-protein': 'v', 'intergenic': '*'}
    facecolor_marker_dict = {'synonymous': white, 'non-synonymous': black, 'non-protein': black, 'intergenic': black}
    facecolor_dict = {True: red, False: blue}
    facecolor_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}
    label_dict = {'synonymous': 'Protein coding - synonymous mutation',
                  'non-synonymous': 'Protein coding - non-synonymous mutation',
                  'non-protein': 'Non-protein coding (rRNA, tRNA etc.)',
                  'intergenic': 'Intergenic (no annotated function)'}

    d['marker'] = sp_df['annotation'].map(marker_dict).values
    d['s'] = -np.log10(
        sp_df[kwargs['pval_col']]).values  # to not divide by max because than it is species subjective
    # d['facecolor'] = sp_df['annotation'].map(facecolor_dict).values
    d['facecolor'] = (sp_df['Coef'] > 0).map(facecolor_dict).values
    d['edgecolor'] = black
    d['linewidths'] = 1
    d['alpha'] = 0.5

    if kwargs['legend_elements'] is None:
        legend_mutation = [Line2D([0], [0], linewidth=0, label=label_dict[k], alpha=d['alpha'],
                                  markersize=10, marker=marker_dict[k], markeredgewidth=d['linewidths'],
                                  markerfacecolor=facecolor_marker_dict[k], markeredgecolor=d['edgecolor']) for k in
                           marker_dict.keys()]

        legend_coefficient = [Line2D([0], [0], linewidth=0, label=facecolor_label_dict[k], alpha=d['alpha'],
                                     markersize=10, marker='s', markeredgewidth=0,
                                     markerfacecolor=facecolor_dict[k], markeredgecolor=white) for k in
                              facecolor_dict.keys()]

        legend = legend_mutation + legend_coefficient
    else:
        legend = kwargs['legend_elements']

    return sp_df, d, legend


def text_func(df: pd.DataFrame, **kwargs):
    text = f"{df.shape[0]:,} SNPs in {df.index.get_level_values('Species').unique().shape[0]:,} species\n" + \
           f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
           f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,} SNPs passed {kwargs['pval_cutoff']:.2e} {kwargs['pval_col'].replace('_', ' ').replace('Pval', 'p-value')} cutoff"
    return text


body_site = 'Oral'

base_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas/PNP3_mwas_{body_site.lower()}_0months_subtraction'
mwas_fname = os.path.join(base_dir, 'mb_gwas.h5')
data_fname = None  # os.path.join(base_dir, 'mb_gwas_data.h5')
annotations_df = None  # pd.read_hdf(os.path.join(base_dir, 'annotations', 'snps_sequences.h5'))[['gene']].rename(columns={'gene': 'text'})
output_dir = os.path.join(base_dir, 'plots')

mwas_plots.run(mwas_fname=mwas_fname, data_fname=data_fname, annotations_df=annotations_df, pval_col='Global_FDR', pval_cutoff=0.1, # input
               out_dir=output_dir, fontsize=10, dpi=200,  # output
               manhattan_draw_func=None, manhattan_text_func=None)
