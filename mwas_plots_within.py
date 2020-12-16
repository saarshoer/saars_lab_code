import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_1_VALUE
from LabData.DataAnalyses.MBSNPs.Plots.manhattan_plot import draw_manhattan_plot


def draw_volcano_plot(df, title='volcano plot', figsize=(10, 6), out_file='volcano_plot',
                      pval_col='pval', pval_cutoff=0.05):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.tick_params(axis='both', which='minor', length=0)
    ax.grid(which='major', axis='x', color='grey', lw=0.4)
    ax.grid(which='major', axis='y', color='grey', lw=0.4)

    # calculations
    coef = df['Coef'].values*MAF_1_VALUE
    pvals = df[pval_col].values

    # volcano
    volcano = ax.scatter(coef, pvals, c=df['N'], alpha=.8, cmap='plasma_r', edgecolor=None, s=10)
    # x
    abs_x_max = max([abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1])])
    ax.set_xlim(-abs_x_max, abs_x_max)
    ax.set_xlabel('coefficient')
    # y
    ax.set_yscale('log')
    ax.set_ylim(1000, pvals.min()/1000)
    ax.set_ylabel('p-value')

    # color bar
    plt.colorbar(mappable=volcano, ax=ax, orientation='vertical', fraction=.05, label='number of samples')

    # significance
    ax.axhline(y=pval_cutoff, linestyle='--', color='red')
    ax.annotate(xy=(-abs_x_max, pval_cutoff), verticalalignment='bottom', s='significance', color='red')

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def draw_qq_plot(df, title='qq plot', figsize=(6, 6), out_file='qq_plot', pval_col='pval'):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.tick_params(axis='both', which='minor', length=0)

    # calculations
    actual_pvals = df['Pval'].dropna().sort_values()  # Pval is intentionally hard coded because it should be raw values
    n = len(actual_pvals)
    expected_pvals = np.arange(n)/n

    # qq
    black = (0.0, 0.0, 0.0, 0.7)
    ax.scatter(x=expected_pvals, y=actual_pvals, marker='o', s=10, facecolors=black, edgecolors=None)
    # actual = expected
    ax.plot(expected_pvals, expected_pvals, '--', color='darkgray')
    # x
    ax.set_xlabel('expected p-value')
    ax.set_xscale('log')
    ax.set_xlim(1, expected_pvals[1])  # 0 fails so has to be the next biggest thing
    # y
    ax.set_ylabel('actual p-value')
    ax.set_yscale('log')
    ax.set_ylim(1000, df[pval_col].min()/1000)

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def run(mwas_fname, annotations_df=None, pval_col='Global_Bonferroni', pval_cutoff=0.05,  # input
        out_dir='.', fontsize=12, dpi=200,  # output
        manhattan_draw_func=None, manhattan_text_func=None):

    rcParams['font.size'] = fontsize
    rcParams['savefig.dpi'] = dpi

    mwas_df = pd.read_hdf(mwas_fname)
    if annotations_df is not None:
        mwas_df = mwas_df.join(annotations_df, on=['Species', 'Contig', 'Position'])
    mwas_df['Pval'] = mwas_df['Pval'].replace(to_replace=0, value=1e-300)  # for qq plot
    mwas_df[pval_col] = mwas_df[pval_col].replace(to_replace=0, value=1e-300)  # for all the rest

    for y, y_df in mwas_df.groupby('Y'):
        y_out_dir = os.path.join(out_dir, y)
        os.makedirs(y_out_dir, mode=0o744, exist_ok=True)

        draw_qq_plot(df=y_df, title=y, out_file=os.path.join(y_out_dir, f'qq_{y}'), pval_col=pval_col)

        draw_volcano_plot(df=y_df, title=y, out_file=os.path.join(y_out_dir, f'volcano_{y}'),
                          pval_col=pval_col, pval_cutoff=pval_cutoff)

        draw_manhattan_plot(df=y_df, title=y, out_file=os.path.join(y_out_dir, f'manhattan_{y}'),
                            draw_func=manhattan_draw_func, text_func=manhattan_text_func,
                            pval_col=pval_col, pval_cutoff=pval_cutoff)


if __name__ == '__main__':

    def mutation_type(annotations_path):
        annotations_df = pd.read_hdf(annotations_path)
        # annotation handling is done with contig without part while mwas data frame is contig with part
        annotations_df = annotations_df.rename(columns={'Contig_with_part': 'Contig'})
        annotations_df = annotations_df.reset_index('Contig', drop=True).set_index('Contig', append=True). \
            reorder_levels(['Species', 'Contig', 'Position'])
        # annotations that can be deduced (so not saved in the hdf files)
        annotations_df['NonSymMutation'] = annotations_df['MajorAA'] != annotations_df['MinorAA']
        annotations_df['annotation'] = annotations_df['NonSymMutation'].map(
            {False: 'synonymous', True: 'non-synonymous'})
        annotations_df.loc[annotations_df['GeneDistance'] < 0, 'annotation'] = 'intergenic'
        annotations_df.loc[(annotations_df['GeneDistance'] >= 0) & (annotations_df['feature'] != 'CDS'), 'annotation'] = \
            'non-protein'  # (df['feature'] != 'CDS') does not cover all options
        # the only thing needed moving forward
        annotations_df = annotations_df['annotation']

        return annotations_df

    def draw_func_annotated(sp_df: pd.DataFrame, **kwargs):
        d = dict()

        black = (0.0, 0.0, 0.0, 0.7)
        white = (1.0, 1.0, 1.0, 0.7)
        blue = (0.0, 0.0, 1.0, 0.7)
        red = (1.0, 0.0, 0.0, 0.7)

        marker_dict = {'synonymous': 'o', 'non-synonymous': 'o', 'non-protein': 'v', 'intergenic': '*'}
        facecolor_marker_dict = {'synonymous': white, 'non-synonymous': black,
                                 'non-protein': black, 'intergenic': black}
        facecolor_dict = {True: red, False: blue}
        facecolor_label_dict = {True: 'positive coefficient', False: 'negative coefficient'}
        label_dict = {'synonymous': 'Protein coding - synonymous mutation',
                      'non-synonymous': 'Protein coding - non-synonymous mutation',
                      'non-protein': 'Non-protein coding (rRNA, tRNA etc.)',
                      'intergenic': 'Intergenic (no annotated function)'}

        d['marker'] = sp_df['annotation'].map(marker_dict).values
        d['s'] = -np.log10(
            sp_df[kwargs['pval_col']]).values  # to not divide by max becaue than it is species subjective
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

    def text_func_annotated_between(df: pd.DataFrame, **kwargs):
        text = f"{df.shape[0]:,} SNPs in {df.index.get_level_values('Species').unique().shape[0]:,} species\n" + \
               f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP\n" + \
               f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,} SNPs passed {kwargs['pval_cutoff']} {kwargs['pval_col'].replace('_', ' ')} cutoff"
        return text

    def text_func_annotated_within(df: pd.DataFrame, **kwargs):
        text = f"{sum(df[kwargs['pval_col']] <= kwargs['pval_cutoff']):,}/{df.shape[0]:,} SNPs passed {kwargs['pval_cutoff']} {kwargs['pval_col'].replace('_', ' ')} cutoff\n" + \
               f"N = between {df['N'].min():,} and {df['N'].max():,} samples per SNP"
        return text

    # def antibiotics():
    run_type = 'within_species'

    input_path = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed', run_type)
    output_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs', run_type)
    jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs/'
    annotations_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/snps_maf_codon_annotations.h5'

    # annotations
    print('starting annotations')
    annotations_df = mutation_type(annotations_path)
    print('finished annotations')

    for file in glob.glob(os.path.join(input_path, 'SGB_*.h5')):  # 9702, 10068
        sgb = 'SGB_' + file.split('SGB_')[1].split('.')[0]
        if len(glob.glob(os.path.join(output_path, sgb, '*'))) != 3:
            try:
                kwargs = {'mwas_fname': file, 'annotations_df': annotations_df, 'out_dir': output_path,
                          'manhattan_draw_func': draw_func_annotated,
                          'manhattan_text_func': text_func_annotated_between if run_type == 'between_species' else text_func_annotated_within}
                run(**kwargs)
            except:
                print(f'{sgb} failed')

    print('done')
