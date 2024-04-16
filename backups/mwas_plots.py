import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import rcParams
from itertools import combinations
from LabQueue.qp import qp, fakeqp
from collections.abc import Iterable
from matplotlib.gridspec import GridSpec
from statannot import add_stat_annotation
from scipy.stats import pearsonr, spearmanr
from LabUtils.addloglevels import sethandlers
from matplotlib.colors import Normalize, LogNorm
from LabData.DataLoaders.MBSNPLoader import MAF_1_VALUE
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
from LabData.DataAnalyses.MBSNPs.Plots.manhattan_plot import draw_manhattan_plot
from LabData.DataAnalyses.MBSNPs.mwas_annots import add_sequences, choose_contig_type

ilegal_chars = [':', '*', '/', '(', ')']

mafft_exe = '/net/mraid20/export/genie/Bin/mafft/mafftdir/bin/mafft'
mafft_binaries = '/net/mraid20/export/genie/Bin/mafft/mafftdir/libexec/'

tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')['Species']


def draw_volcano_plot(df, title='volcano plot', figsize=(10, 6), out_file='volcano_plot',
                      coef_col='Coef', pval_col='Pval', pval_cutoff=0.05):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.tick_params(axis='both', which='minor', length=0)
    ax.grid(which='major', axis='x', color='grey', lw=0.4)
    ax.grid(which='major', axis='y', color='grey', lw=0.4)

    # calculations
    coef = df[coef_col].values*MAF_1_VALUE
    pvals = df[pval_col].values

    # volcano
    volcano = ax.scatter(coef, pvals, c=df['N'], alpha=.8, cmap='plasma_r', edgecolor=None, marker='o', s=10)
    # x
    abs_x_max = max([abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1])])
    ax.set_xlim(-abs_x_max, abs_x_max)
    ax.set_xlabel('coefficient')
    # y
    ax.set_yscale('log')
    min_pval = df[pval_col].min()
    min_pval_in_fig = min_pval * 10 ** (0.1 * np.log10(min_pval))  # Add 10% white space above
    if min_pval_in_fig == 0:  # Overflow
        min_pval_in_fig = min_pval
    ax.set_ylim(1, min_pval_in_fig)
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


def draw_qq_plot(df, title='qq plot', figsize=(6, 6), out_file='qq_plot', pval_col='Pval'):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.tick_params(axis='both', which='minor', length=0)

    # calculations
    df = df.reset_index().drop_duplicates(['Species', 'Contig', 'Position'])
    # because rows can be multiplied if they have multiple genes
    actual_pvals = df[pval_col].dropna().sort_values()  # should be raw values
    n = len(actual_pvals)
    expected_pvals = np.arange(1, n+1)/n

    # qq
    black = (0.0, 0.0, 0.0, 0.7)
    min_pval = min(expected_pvals[0], actual_pvals.iloc[0])
    min_pval_in_fig = min_pval * 10 ** (0.1 * np.log10(min_pval))  # Add 10% white space above
    if min_pval_in_fig == 0:  # Overflow
        min_pval_in_fig = min_pval
    ax.plot([min_pval_in_fig, 1], [min_pval_in_fig, 1], '--', color='darkgray')
    ax.scatter(x=expected_pvals, y=actual_pvals, marker='o', s=10, facecolors=black, edgecolors=None)

    # x
    ax.set_xlabel('expected p-value')
    ax.set_xscale('log')
    ax.set_xlim(1, min_pval_in_fig)
    # y
    ax.set_ylabel('actual p-value')
    ax.set_yscale('log')
    ax.set_ylim(1, min_pval_in_fig)

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def draw_box_plot(df, y_df=None, title='box plot', figsize=(12, 6), out_file='box_plot', maf_col='MAF', per_bin=True):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=(figsize[0] if per_bin else figsize[0]/3, figsize[1]))

    # add info for people without the snp
    if y_df is not None:
        y_col = df.index.get_level_values('Y').unique().values[0]
        y_df = y_df.loc[y_df.index.difference(df.index.get_level_values('SampleName')), y_col] \
            .to_frame().rename(columns={y_col: 'y'})
        df = pd.concat([df.reset_index().set_index('SampleName'), y_df]).reset_index().set_index(df.index.names)

    # calculations
    labels = None if per_bin else ['minor', 'major']
    bins = list(np.arange(0.0, 1.1, 0.1)) if per_bin else [0.0, 0.5, 1.0]
    if y_df is not None:
        bins = bins + [1.1]
        if labels is not None:
            labels = labels + ['no snp']
    df['bin'] = pd.cut(df[maf_col].fillna(1.1*MAF_1_VALUE)/MAF_1_VALUE, bins=bins, labels=labels, include_lowest=True)

    n_colors = 51  # should be an odd number so to have white at zero
    cmap = sns.color_palette(palette='coolwarm', n_colors=n_colors)
    color_indices = (df.groupby('bin')['y'].median().fillna(0) /
                     df.groupby('bin')['y'].median().abs().max() + 1) * 0.5 * (n_colors - 1)
    colors = [cmap[i] for i in color_indices.astype(int)]

    # box plot
    sns.boxplot(x='bin', y='y', palette=colors, data=df, ax=ax)

    xlabels = ax.get_xticklabels()
    d = df['bin'].value_counts()
    d.index = d.index.astype(str)
    d = d.to_dict()
    for l in xlabels:
        l.set_text(f'{l.get_text()}\nn={d[l.get_text()]}')
    ax.set_xticklabels(labels=xlabels, rotation=45)
    ax.set_ylabel(f"change in {df.index.get_level_values('Y')[0]}")

    if per_bin:
        xlabels = ax.get_xticklabels()
        xlabels[0].set_text(xlabels[0].get_text().replace('-0.001', '0.0'))
        xlabels[-1].set_text(xlabels[-1].get_text().replace('(1.0, 1.1]', 'no snp'))
        ax.set_xlabel('Major Allele Frequency')

        r_p, p_p = pearsonr(df[maf_col].dropna()/MAF_1_VALUE, df.dropna(subset=[maf_col])['y'])
        r_s, p_s = spearmanr(df[maf_col].dropna()/MAF_1_VALUE, df.dropna(subset=[maf_col])['y'])
        text = f"pearson: r={r_p:.2f} p={p_p:.2e}\nspearman: r={r_s:.2f} p={p_s:.2e}"
        ax.text(x=0.7, y=0.875, s=text, transform=ax.transAxes)
    else:
        ax.set_xlabel('Allele')

        ax, test_results = add_stat_annotation(ax, x='bin', y='y', data=df, box_pairs=list(combinations(labels, 2)),
                                               test='Mann-Whitney', text_format='star')

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def draw_scatter_plot(df, title='scatter plot', figsize=(12, 6), out_file='scatter_plot', maf_col='MAF'):

    # figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # calculations
    r_p, p_p = pearsonr(df[maf_col]/MAF_1_VALUE, df['y'])
    r_s, p_s = spearmanr(df[maf_col]/MAF_1_VALUE, df['y'])
    corr_text = f"pearson: r={r_p:.2f} p={p_p:.2e}\nspearman: r={r_s:.2f} p={p_s:.2e}"

    # box plot
    sns.scatterplot(x=maf_col, y='y', hue='bin', palette='coolwarm', data=df, ax=ax)
    ax.set_xlabel('Major Allele Frequency')
    ax.set_ylabel(f"delta change in {df.index.get_level_values('Y')[0]}")
    # ax.set_ylabel(f"{df.index.get_level_values('Y')[0]}")
    ax.text(x=0.7, y=0.875, s=corr_text, transform=ax.transAxes)

    # finish
    ax.set_title(title)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


def draw_msa_plot(df, title='msa plot', figsize=(10, 6), out_file='msa_plot',
                  pval_col='Pval', coef_col='Coef', rc=None,
                  mafft_arguments='--localpair  --maxiterate 16 --reorder',
                  nuc=False, drop_cols_without_snp=False):

    # TODO: multiple references with excel
    # TODO: excel color with partial positions

    global rcParams
    rc_backup = rcParams.copy()
    if rc is not None:
        rcParams.update(rc)

    # whether plotting nucleutides or amino acids
    seq_col = 'nuc_seq' if nuc else 'aa_seq'
    major_col = 'MajorAllele' if nuc else 'MajorAA'
    minor_col = 'MinorAllele' if nuc else 'MinorAA'

    # create folder for output files
    out_dir = os.path.dirname(out_file)
    os.makedirs(out_dir, exist_ok=True)

    df['GeneID'] = df['GeneID'].astype(int)
    df[coef_col] = df[coef_col].astype(float)*MAF_1_VALUE
    df[['MajorAllele', 'MinorAllele']] = df[['MajorAllele', 'MinorAllele']].replace({0: 'a', 1: 'c', 2: 'g', 3: 't'})

    seq_file = os.path.join(out_dir, f'{title}_{"nuc" if nuc else "aa"}_sequences.fasta')
    align_file = os.path.join(out_dir, f'{title}_{"nuc" if nuc else "aa"}_sequences_alignment.fasta')
    ref_file = os.path.join(out_dir, f'{title}_{"nuc" if nuc else "aa"}_sequences_reference.fasta')

    # retrieve reference sequences
    if not os.path.exists(seq_file):
        df = add_sequences(choose_contig_type(df, 'Contig_without_part'))  # extremely slow in debug
        # write to file
        df['4fasta'] = '>' + df['GeneID'].astype(str)
        df[['4fasta', seq_col]].drop_duplicates().set_index('4fasta').to_csv(
            os.path.join(out_dir, seq_file), header=False, sep='\n')
        df = df.drop('4fasta', axis=1)

    # multiple sequence alignment
    if not os.path.exists(align_file):
        os.environ['MAFFT_BINARIES'] = mafft_binaries
        if os.path.exists(ref_file):
            cmd = f'echo "" >> {seq_file}'
            os.system(cmd)
            cmd = f'cat {ref_file} >> {seq_file}'
            os.system(cmd)
        cmd = f'{mafft_exe}  {mafft_arguments} {seq_file} > {align_file}'
        os.system(cmd)
    # read from file
    align_order = []
    references_text = dict()
    seqs = SeqIO.parse(open(align_file), 'fasta')
    with open(align_file) as file:
        for seq in seqs:
            try:
                df.loc[df['GeneID'] == int(seq.id), f'aligned_{seq_col}'] = str(seq.seq)
                align_order.append(int(seq.id))
            except ValueError:
                references_text[str(seq.description)] = str(seq.seq)

    # data frames for plots
    text = df.dropna(subset=['GeneID']).drop_duplicates(['GeneID']).reset_index().set_index('GeneID')
    text = text[f'aligned_{seq_col}'].apply(lambda s: pd.Series([s[i] for i in np.arange(len(s))]))
    references_text = pd.DataFrame.from_dict(references_text, orient='index', columns=[f'aligned_{seq_col}'])
    references_text = references_text[f'aligned_{seq_col}'].apply(lambda s: pd.Series([s[i] for i in np.arange(len(s))]))

    # add snp information to adjusted alignment position
    pvals = pd.DataFrame(index=text.index, columns=text.columns)
    coefs = pd.DataFrame(index=text.index, columns=text.columns)

    for snp in df.index:
        gene_id = df.loc[snp, 'GeneID']
        pos_in_seq = int(df.loc[snp, 'GeneDistance']) if nuc else int(np.floor(df.loc[snp, 'GeneDistance']/3))
        seq2align_pos = {text.loc[gene_id, :pos_in_align].replace('-', np.nan).dropna().shape[0]-1: pos_in_align
                         for pos_in_align in text.columns}

        text.loc[gene_id, seq2align_pos[pos_in_seq]] = f'{df.loc[snp, major_col]}/{df.loc[snp, minor_col]}'
        pvals.loc[gene_id, seq2align_pos[pos_in_seq]] = df.loc[snp, pval_col]
        coefs.loc[gene_id, seq2align_pos[pos_in_seq]] = df.loc[snp, coef_col]

    pvals = pvals.loc[align_order] if not drop_cols_without_snp else pvals.loc[align_order].dropna(how='all', axis=1)
    coefs = coefs.loc[align_order, pvals.columns]
    text = text.loc[align_order, pvals.columns]

    references_excel = pd.DataFrame()
    if len(references_text) > 0:

        # xlabels - position in reference gene instead of gene alignment
        align2ref_pos = {pos_in_align: references_text.loc[references_text.index[-1], :pos_in_align]
                         .replace('-', np.nan).dropna().shape[0] for pos_in_align in references_text.columns}
        references_text = references_text.loc[:, pvals.columns]
        references_text.columns = [align2ref_pos[pos_in_align] for pos_in_align in references_text.columns]

        # adjust reference according to excel
        if os.path.exists(ref_file.replace('.fasta', '.xlsx')):
            references_excel = pd.read_excel(ref_file.replace('.fasta', '.xlsx'), header=1)
            references_excel = references_excel.set_index('position').T
            missing_positions = references_text.columns.difference(references_excel.columns)
            if len(missing_positions) > 0:
                references_excel.loc[:, missing_positions] = np.nan
            references_excel = references_excel.loc[:, references_text.columns]
            if 'text' in references_excel.index:
                for pos, s in references_excel.loc['text'].dropna().iteritems():
                    pos_index = references_text.columns.get_loc(pos)
                    pos_index = pos_index[0] if isinstance(pos_index, Iterable) is int else pos_index
                    references_text.iloc[-1, pos_index] = s

    # xlabels - 1 based position numbering
    pvals.columns = coefs.columns = text.columns = [col + 1 for col in pvals.columns]

    # ylabels
    pvals.index = pvals.index.map((df.reset_index().set_index('GeneID')['Species'] +
                                   ' (' + df['GeneID'].astype(int).astype(str).values + ')').to_dict())
    coefs.index = pvals.index.copy()

    conservation = text.apply(lambda col: int(100*(col == col.value_counts().index[0]).sum()/col.shape[0])). \
        reset_index().rename(columns={'index': 'Position in Gene', 0: 'Allele\nconservation'})

    # plot
    for col in [pval_col, coef_col]:

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.axis('off')
        n_con = 3
        n_seq = pvals.shape[0]
        n_ref = max(1, references_text.shape[0])
        n_pos = pvals.shape[1]
        gs = GridSpec(n_con+n_seq+n_ref+1, n_pos+1, hspace=1)
        ax_conservation: plt.Axes = fig.add_subplot(gs[:n_con, :n_pos])
        ax_heatmap: plt.Axes = fig.add_subplot(gs[n_con:n_con+n_seq, :n_pos])
        ax_cbar: plt.Axes = fig.add_subplot(gs[n_con:n_con+n_seq, n_pos])
        ax_reference: plt.Axes = fig.add_subplot(gs[n_con+n_seq:n_con+n_seq+n_ref, :n_pos])

        # ax_conservation
        sns.barplot(x='Position in Gene', y='Allele\nconservation', data=conservation,
                    color='lightgrey',
                    ax=ax_conservation)
        ax_conservation.set_title(title)
        ax_conservation.tick_params(top=False, bottom=False, left=False, right=True,
                                    labeltop=False, labelbottom=False, labelleft=False, labelright=True)
        ax_conservation.set_xlabel('')
        ax_conservation.yaxis.set_label_position('right')
        ax_conservation.set_ylim(0, 100)

        # ax_heatmap
        if col == pval_col:
            data = pvals.fillna(1)
            norm = LogNorm(max(1e-300, pvals.min().min()), 1, clip=False)  # smaller number ends up to be zero
            cmap = 'afmhot'
            label = 'p-value'
        else:
            data = coefs.fillna(0)
            norm = Normalize(-coefs.abs().max().max(), coefs.abs().max().max(), clip=False)
            cmap = 'bwr'
            label = 'Coefficient'

        sns.heatmap(data=data, norm=norm, cmap=cmap,
                    cbar=True, cbar_kws={'label': label}, cbar_ax=ax_cbar,
                    annot=text, fmt='',
                    ax=ax_heatmap)
        ax_heatmap.tick_params(top=False, bottom=False, left=False, right=False,
                               labeltop=False, labelbottom=False if len(references_text) > 0 else True,
                               labelleft=True, labelright=False)
        ax_heatmap.set_yticks(np.arange(data.shape[0]) + 0.5)
        ax_heatmap.set_yticklabels(data.index, rotation=0)
        ax_heatmap.set_ylabel('')

        # ax_reference
        if len(references_text) > 0:

            if len(references_excel) > 0:
                label2int = {l: i for i, l in enumerate(references_excel.loc['label'].fillna('').unique())}
                data = references_excel.loc[['label']].fillna('').replace(label2int)
                cmap = references_excel.loc['color'].fillna('grey').tolist() if 'color' in references_excel.index else \
                       sns.color_palette('Pastel1', len(label2int))
                xticklabels = list(references_excel.loc['label'].fillna('') +
                                   references_excel.columns.astype(str).str.pad(4))
            else:
                data = references_text.isna()
                cmap = ['grey']
                xticklabels = references_text.columns

            sns.heatmap(data=data,
                        cbar=False, cmap=cmap,
                        annot=references_text, fmt='',
                        ax=ax_reference)
            ax_reference.tick_params(top=False, bottom=False, left=False, right=False,
                                     labeltop=False, labelbottom=True, labelleft=True, labelright=False)
            ax_reference.set_xlabel('Position in Reference Gene')
            ax_reference.set_xticks(np.arange(len(xticklabels)) + 0.5)
            ax_reference.set_xticklabels(xticklabels, rotation=90)
            ax_reference.set_yticks(np.arange(len(references_text.index)) + 0.5)
            ax_reference.set_yticklabels(references_text.index.str.split(' #').str[0], rotation=0)
        else:
            ax_reference.axis('off')
            ax_heatmap.set_xlabel('Position in Gene Alignment')
            ax_heatmap.set_xticks(np.arange(data.shape[1]) + 0.5)
            ax_heatmap.set_xticklabels(data.columns, rotation=90)
        # finish
        plt.savefig(f'{out_file}_{"nuc" if nuc else "aa"}_{col}', bbox_inches='tight')
        plt.close()

    rcParams.clear()
    rcParams.update(rc_backup)


def draw_msa_plot_from_jupyter(d, jobs_dir):
    # for unexplainable reason, mafft is unable to run when draw_msa_plot is imported to jupyter
    # (it is able to run directly from jupyter)
    # so this is a work around
    os.makedirs(jobs_dir, exist_ok=True)
    os.chdir(jobs_dir)
    try:  # only needs to be done once
        sethandlers(file_dir=jobs_dir)
    except:
        pass

    with qp(jobname='msa', _tryrerun=False, _suppress_handlers_warning=True) as q:
        q.startpermanentrun()
        tkttores = {}

        for job_name, kwargs in d.items():
            kwargs['rc']['font.size'] = 10
            tkttores[job_name] = q.method(draw_msa_plot, kwargs=kwargs, _job_name=job_name)

        for k, v in tkttores.items():
            q.waitforresult(v, _assert_on_errors=False)


def run(mwas_fname=None, data_fname=None, annotations_df=None, y_df=None,  # input
        maf_col='MAF', coef_col='Coef', pval_col='Global_Bonferroni', pval_cutoff=0.05,  pval_func=None,  # input
        out_dir='.', fontsize=10, dpi=200,  # output
        manhattan_draw_func=None, manhattan_text_func=None):

    rcParams['font.size'] = fontsize
    rcParams['savefig.dpi'] = dpi

    if data_fname is not None:
        if type(data_fname) is str:
            data_df = pd.read_hdf(data_fname)
        else:
            data_df = data_fname
        data_df['Position'] = data_df.index.get_level_values('Position').astype(int)
        data_df = data_df.reset_index('Position', drop=True).set_index('Position', append=True)
        if annotations_df is not None:
            data_df = data_df.join(annotations_df, on=['Y', 'Species', 'Contig', 'Position'])#.dropna(subset=[pval_col])
            # multiply rows in caller if there are multiple matches in other

        data_out_dir = os.path.join(out_dir, 'data')
        os.makedirs(data_out_dir, mode=0o744, exist_ok=True)

        for (y, species, contig, position), snp_df in data_df.groupby(['Y', 'Species', 'Contig', 'Position']):
            y_legal = y
            for char in ilegal_chars:
                y_legal = y_legal.replace(char, '')
            title = f"{tax_df.loc[species].split('s__')[-1]} and {y}\n{species} {contig} {position}"

            snp_df = snp_df.rename(columns={position: maf_col})

            draw_box_plot(snp_df.copy(), y_df=y_df.copy() if y_df is not None else None, title=title,
                          out_file=os.path.join(data_out_dir, f'box_plot_{species}_{contig}_{position}_{y_legal}'),
                          maf_col=maf_col, per_bin=False)

            # draw_scatter_plot(snp_df.copy(), title=title,
            #                   out_file=os.path.join(data_out_dir, f'scatter_plot_{species}_{contig}_{position}_{y_legal}'))

    if mwas_fname is not None:
        if type(mwas_fname) is str:
            mwas_df = pd.read_hdf(mwas_fname)
        else:
            mwas_df = mwas_fname
        if annotations_df is not None:
            mwas_df = mwas_df.join(annotations_df, on=['Y', 'Species', 'Contig', 'Position'])#.dropna(subset=[pval_col])
            # multiply rows in caller if there are multiple matches in other

        # smallest round non-zero
        min_value = mwas_df.loc[mwas_df[pval_col] != 0, pval_col].min()
        min_value = min_value if min_value > 1e-300 else 1e-300  # smaller number ends up to be zero
        mwas_df[pval_col] = mwas_df[pval_col].clip(lower=10**-(math.ceil(-np.log10(min_value)/10)*10))
        # if there is co-linearity between features the p-value will be missing
        condition = (mwas_df[pval_col].isna()) & (mwas_df[coef_col] == 0)
        if sum(condition) > 0:
            mwas_df.loc[condition, pval_col] = 1

        for y, y_df in mwas_df.groupby('Y'):
            y_legal = y
            for char in ilegal_chars:
                y_legal = y_legal.replace(char, '')
            y_out_dir = os.path.join(out_dir, y_legal)
            ys_pval_cutoff = pval_cutoff if pval_func is None else pval_func(pval_col, pval_cutoff, y_df)
            os.makedirs(y_out_dir, mode=0o744, exist_ok=True)

            title = f"{y}\n{tax_df.loc[y].split('s__')[-1]}" if 'SGB' in y or 'Rep' in y else y

            draw_qq_plot(df=y_df, title=title, out_file=os.path.join(y_out_dir, f'qq_{y_legal}'), pval_col=pval_col)

            draw_volcano_plot(df=y_df, title=title, out_file=os.path.join(y_out_dir, f'volcano_{y_legal}'),
                              coef_col=coef_col, pval_col=pval_col, pval_cutoff=ys_pval_cutoff)

            for col in [pval_col, coef_col]:
                draw_manhattan_plot(df=y_df, title=title,
                                    out_file=os.path.join(y_out_dir, f'manhattan_{y_legal}_{col}'),
                                    draw_pval=col==pval_col,
                                    draw_func=manhattan_draw_func, text_func=manhattan_text_func,
                                    coef_col=coef_col, pval_col=pval_col, pval_cutoff=ys_pval_cutoff, rc=rcParams)


if __name__ == '__main__':
    base_path = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/'
    run_dir = '/net/mraid20/export/genie/LabData/Analyses/saarsh/anti_mwas_within'

    annots = pd.read_hdf(os.path.join(run_dir, 'snps_gene_annotations_short.h5'))
    gene = 'ssrA'

    draw_msa_plot(**{'df': annots[(annots['text'] == gene)],
                                  'pval_col': 'Pval',
                                  'coef_col': 'Coef',
                                  'nuc': True,
                                  'drop_cols_without_snp': True,
                                  'title': gene,
                                  'out_file': os.path.join(base_path, 'figs', f'within_genes', gene,
                                                           f'debug_{gene}'),
                                  'rc': rcParams})
