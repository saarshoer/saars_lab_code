import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from LabData.DataLoaders.MBSNPLoader import MAF_1_VALUE
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
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
    volcano = ax.scatter(coef, pvals, c=df['N'], alpha=.8, cmap='plasma_r', edgecolor=None, marker='o', s=10)
    # x
    abs_x_max = max([abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1])])
    ax.set_xlim(-abs_x_max, abs_x_max)
    ax.set_xlabel('coefficient')
    # y
    ax.set_yscale('log')
    ax.set_ylim(1, pvals.min())
    ax.set_ylabel('p-value')

    # text
    if 'text' in df.columns:
        for (x, y, t) in zip(coef, pvals, df['text']):
            if type(t) is str:
                ax.text(x, y, t)

    # color bar
    plt.colorbar(mappable=volcano, ax=ax, orientation='vertical', fraction=.05, label='number of samples')

    # significance
    ax.axhline(y=pval_cutoff, linestyle='--', color='red')
    # ax.annotate(xy=(-abs_x_max, pval_cutoff), xycoords='data', text='significance', color='red')

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def draw_qq_plot(df, title='qq plot', figsize=(6, 6), out_file='qq_plot', pval_col='pval'):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.tick_params(axis='both', which='minor', length=0)

    # calculations
    df = df.reset_index().drop_duplicates(['Species', 'Contig', 'Position'])
    # because rows can be multiplied if they have multiple genes
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
    x_max = expected_pvals[1] if len(expected_pvals) > 1 else 10**-1
    ax.set_xlim(1, x_max)  # 0 fails so it has to be the next smallest thing, but sometimes 0 is all you have
    # y
    ax.set_ylabel('actual p-value')
    ax.set_yscale('log')
    # intentionally pval_col and not 'Pval' so to be the same as in other graphs
    ax.set_ylim(1, df[pval_col].min())

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def run(mwas_fname, annotations_df=None, pval_col='Global_Bonferroni', pval_cutoff=0.05,  # input
        out_dir='.', fontsize=10, dpi=200,  # output
        manhattan_draw_func=None, manhattan_text_func=None):

    rcParams['font.size'] = fontsize
    rcParams['savefig.dpi'] = dpi

    mwas_df = pd.read_hdf(mwas_fname)
    if annotations_df is not None:
        mwas_df = mwas_df.join(annotations_df, on=['Species', 'Contig', 'Position']).dropna(subset=[pval_col])
        # multiply rows in caller if there are multiple matches in other
    # for qq plot - smallest round non-zero
    min_value = mwas_df.loc[mwas_df['Pval'] != 0, 'Pval'].min()
    min_value = min_value if min_value > 1e-319 else 1e-319  # smaller number ends up to be zero
    mwas_df['Pval'] = mwas_df['Pval'].replace(to_replace=0, value=10**-(math.ceil(-np.log10(min_value)/10)*10))
    # for all the rest - smallest round non-zero
    min_value = mwas_df.loc[mwas_df['Pval'] != 0, 'Pval'].min()
    min_value = min_value if min_value > 1e-319 else 1e-319  # smaller number ends up to be zero
    mwas_df[pval_col] = mwas_df[pval_col].replace(to_replace=0, value=10**-(math.ceil(-np.log10(min_value)/10)*10))

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')['Species']

    for y, y_df in mwas_df.groupby('Y'):
        y_out_dir = os.path.join(out_dir, y)
        os.makedirs(y_out_dir, mode=0o744, exist_ok=True)

        title = f"{y}\n{tax_df.loc[y].split('s__')[-1]}" if 'SGB' in y else y

        draw_qq_plot(df=y_df, title=title, out_file=os.path.join(y_out_dir, f'qq_{y}'), pval_col=pval_col)

        draw_volcano_plot(df=y_df, title=title, out_file=os.path.join(y_out_dir, f'volcano_{y}'),
                          pval_col=pval_col, pval_cutoff=pval_cutoff)

        draw_manhattan_plot(df=y_df, title=title, out_file=os.path.join(y_out_dir, f'manhattan_{y}'),
                            draw_func=manhattan_draw_func, text_func=manhattan_text_func,
                            pval_col=pval_col, pval_cutoff=pval_cutoff, rc=rcParams)
