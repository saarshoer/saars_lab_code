import os
import glob
import numpy as np
import pandas as pd

from numpy.linalg import LinAlgError
from statsmodels.api import Logit, OLS
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from statsmodels.tools.sm_exceptions import PerfectSeparationError

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import rcParams, cm, colors

from analysis import segal_name
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE, MAF_1_VALUE

rcParams['font.size'] = 12
rcParams['figure.figsize'] = (rcParams['figure.figsize'][0] * 1.5, rcParams['figure.figsize'][1] * 1.5)  # figure size in inches
rcParams['savefig.dpi'] = 200  # figure dots per inch or 'figure'
rcParams['savefig.format'] = 'pdf'  # {png, ps, pdf, svg}
rcParams['savefig.bbox'] = 'tight'

assemblies = pd.read_csv('/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/Unicorn/URS/URS_Build/full_metadata_filtered_w_clusts_reps.csv', index_col=0)
assemblies.loc[assemblies['Source'] == 'Segal', 'SampleName'] = assemblies[assemblies['Source'] == 'Segal'].index.str.split('/').str[-3].str.split('-').str[-1]
assemblies.loc[assemblies['Source'] == 'Segata', 'SampleName'] = assemblies[assemblies['Source'] == 'Segata'].index.str.split('/').str[-1].str.split('__bin').str[0].str.replace('.', '')
assemblies.loc[~assemblies['Source'].isin(['Segal', 'Segata']), 'SampleName'] = assemblies[~assemblies['Source'].isin(['Segal', 'Segata'])].index.str.split('/').str[-1].str.split('#').str[0].str.split('_ASM').str[0].str.replace('.fa', '').str.replace('.fna', '').str.replace('.', '_')
assemblies['SampleName'] = assemblies['SampleName'].str.split('_v').str[0].str.replace('sample_', '')
assemblies = assemblies.reset_index()

# snps
maf = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/antibiotics/{}/pcs_covariate/mb_snp_maf_{}.h5'
strain_variability = '/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/mb_sns_strain_variability_R3.h5'

# pangenomes
dists_mat16 = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/Unicorn/URS/URS_Build/dists/between_person/dists_mat16.dat'
gene_presence_absence = '/net/mraid20/export/genie/LabData/Data/Annotations/Segal_annots_all_assemblies/{}/roary_3.13.0/gene_presence_absence.csv'

mash_dist = None#np.random.uniform(low=0, high=1, size=(241118, 241118))


def get_distance(species, metric, maximal_imputation, df=None):
    if metric == 'diss':

        # get data
        dist = pd.read_hdf(os.path.join('data_frames', 'mb_dists.h5'), key=species)
        dist = dist.xs('baseline', level='time_point1').xs('baseline', level='time_point2')
        dist = dist.reset_index().pivot_table(index='SampleName2', columns='SampleName1', values='dissimilarity')

        # limit to samples that have at most maximal_filling percentage of existing dissimilarity measurements
        samples_mask = dist.columns[dist.isna().sum() < dist.shape[0] * maximal_imputation].values
        dist = dist.loc[samples_mask, samples_mask]

        # fill missing values
        mask = dist.isna()
        dist = dist.apply(lambda row: row.fillna(row.median()))

    elif metric == 'mash':

        # get data
        global mash_dist
        if mash_dist is None:
            mem = np.memmap(dists_mat16, dtype=np.float16, mode='r', shape=(len(assemblies), len(assemblies)))
            mash_dist = np.array(mem)
            del mem
        genomes = assemblies.loc[assemblies['cluster_s'] == int(species.split('_')[-1]), 'SampleName']
        dist = pd.DataFrame(mash_dist[genomes.index][:, genomes.index], index=genomes.values, columns=genomes.values)
        mask = None

    else:
        mask = df.isna()
        df = df.fillna(1)  # only happens for maf
        dist = pd.DataFrame(squareform(pdist(df, metric=metric)), index=df.index, columns=df.index)

    return dist, mask


def get_linkage(dist, method):
    return linkage(squareform(dist, force='tovector', checks=False), method=method, optimal_ordering=True)


def find_strains(dist, linkage_, relative_change_threshold, min_samples_per_strain, max_samples_per_strain):

    colors_ = [mpl.colors.rgb2hex(c) for c in cm.Set1.colors[:-1]]
    color_names = ['re', 'bl', 'gr', 'pu', 'or', 'ye', 'br', 'pi']
    color_names = [f'{c}3' for c in color_names] + [f'{c}2' for c in color_names] + color_names
    color_names = dict(zip(colors_ + colors_ + colors_, color_names))
    color_names['black'] = 'black'

    j = 1
    link_colors = {}

    for i, link in reversed(list(enumerate(linkage_))):
        daddy = np.where(linkage_[:, :2] == (i + dist.shape[0]))[0]
        daddy = daddy[0] if daddy.shape[0] == 1 else linkage_.shape[0] - 1

        relative_change = 1 - link[2] / linkage_[daddy, 2]

        if i != daddy and link_colors[daddy + dist.shape[0]] != 'black':
            link_colors[i + dist.shape[0]] = link_colors[daddy + dist.shape[0]]
        elif (link[3] < dist.shape[0] * max_samples_per_strain) and \
             (link[3] >= min_samples_per_strain) and\
             (relative_change >= relative_change_threshold):
            link_colors[i + dist.shape[0]] = colors_[j % len(colors_)]
            j = j + 1
        else:
            link_colors[i + dist.shape[0]] = 'black'

    return link_colors, color_names


def get_data(species, snps):
    if snps:
        file = maf.format(study, species)
        data = pd.read_hdf(file).replace(MAF_MISSING_VALUE, np.nan) / MAF_1_VALUE
        data = data.T

    else:
        file = gene_presence_absence.format(species)

        if species == 'Rep_449':
            all_chunks = []
            for chunk in pd.read_csv(file, dtype=str, index_col=0, chunksize=10000):
                all_chunks.append((~chunk.iloc[:, 13:].isna()))
            del chunk
            data = pd.concat(all_chunks)
            del all_chunks
        else:
            data = (~pd.read_csv(file, dtype=str, index_col=0).iloc[:, 13:].isna())

        data.columns = data.columns.str.split('_v').str[0].str.replace('sample_', '')
        freq = data.sum(axis=1) / data.shape[1] * 100
        data = data[(freq > 10) & (freq <= 90)]  # shell genes

    return data


def get_statistics(a, covariates, df, sample_colors, min_samples_per_strain, flip=False):

    def test(c1, c2=None):

        if c1 == 'continuous':
            model = OLS
            data4model = sample_colors_df.copy()
        elif len(df[a].unique()) == 2:
            model = Logit
            data4model = sample_colors_df == c1 if c2 is None else \
                sample_colors_df[sample_colors_df != 'black'] == c1
        else:
            model = OLS
            data4model = sample_colors_df == c1 if c2 is None else \
                sample_colors_df.loc[sample_colors_df.isin([c1, c2])] == c1

        data4model = data4model.to_frame('strain').join(
            df[list(set(covariates + [a]))], how='inner').assign(constant=1)
        data4model = data4model.replace({False: 0, True: 1, 'Female': 0, 'Male': 1, 'Unknown': np.nan})

        if flip:
            x = a
            y = 'strain'
        else:
            x = 'strain'
            y = a

        try:
            results = model(data4model[y], data4model.loc[:, data4model.columns != y],
                            missing='drop').fit(disp=False)
        except LinAlgError:
            results = model(data4model[y], data4model.loc[:, data4model.columns != y],
                            missing='drop').fit(method='nm', disp=False)
        except PerfectSeparationError:
            return 1, 0
        r2 = results.prsquared if (c1 != 'continuous') & (len(df[a].unique()) == 2) else results.rsquared

        if r2 > 0:
            return results.params[x], results.pvalues[x]
        else:
            return 0, 1

    sample_colors_df = pd.Series([v[1] for v in sample_colors.values()], index=[v[0] for v in sample_colors.values()])
    sample_colors_vc = sample_colors_df.value_counts()
    if 'black' in sample_colors_vc.index:
        n_black = sample_colors_vc.loc['black']
        sample_colors_vc = sample_colors_vc[sample_colors_vc.index != 'black']
    else:
        n_black = 0
    n_colors = len(sample_colors_vc)

    if n_colors == 0 or (n_colors == 1 and n_black < min_samples_per_strain):
        # nothing
        p = 1
        l = a

    elif n_colors == 2 and n_black < min_samples_per_strain:
        # 1 vs 1
        color1 = sample_colors_vc.index[0]
        color2 = sample_colors_vc.index[1]
        r, p = test(color1, color2)
        l = f'{a} ~ {color1} c={r:.2f} p={p:.0e}'

    else:
        # 1 vs all
        all_p = []
        anno_label = f'{a} ~'
        for color in (sample_colors_vc.index if sample_colors_vc.shape[0] < 20 else ['continuous']):
            r, p = test(color)
            all_p.append(p)
            anno_label = anno_label + (f' {color} c={r:.2f} p={p:.0e}' if n_colors > 1 else f' {color} c={r:.2f} p={p:.0e}')
        p = min(all_p)
        l = anno_label

    return p, l


def filter_data(df, sample_colors, min_samples_per_strain, suffix):
    results = pd.DataFrame(index=df.index, columns=['pvalue', 'label'])
    df = df.T
    for anno in df.columns:  # genes
        results.loc[anno, ['pvalue', 'label']] = get_statistics(anno, [], df, sample_colors, min_samples_per_strain, True)
    results['p_Bonfferoni'] = (results['pvalue'] * results.shape[0]).clip(upper=1)
    df = df[results[results['p_Bonfferoni'] < 0.05].index]
    df = df.T
    mask = df.isna()

    # saving
    subdir = 'snps' if snps else 'genes'
    if not os.path.exists(os.path.join('figs', subdir, 'data')):
        os.makedirs(os.path.join('figs', subdir, 'data'))
    results.to_csv(os.path.join('figs', subdir, 'data', f'{species}{suffix}.csv'))

    return df, mask


def get_colors(df, cmaps):

    known_colors = {True: 'black',  False: 'white',
                    'Male': 'lightskyblue', 'Female': 'lightpink', 'Unknown': 'black'}

    # replace values with corresponding colors
    for a in df.columns:
        if len(set(df[a].dropna()) - set(known_colors.keys())) > 0:  # most probably continuous
            m = cm.ScalarMappable(
                norm=colors.Normalize(vmin=df[a].min(), vmax=df[a].max()),
                cmap=cmaps[a] if a in cmaps.keys() else None)
            df[a] = df[a].apply(m.to_rgba).to_frame()
        else:  # discrete
            df[a] = df[a].replace(known_colors)

    return df


def plot(df, mask, distance, snps, significant, sig_genes,
         row_linkage, col_linkage, link_colors,
         row_colors, col_colors):

    if distance:
        cmap = 'autumn'
    else:
        if snps:
            cmap = 'Blues'
            rep = None
        else:
            cmap = colors.ListedColormap(['gold', 'white'])
            rep = assemblies.loc[(assemblies['cluster_s'] == int(species.split('_')[-1])) & assemblies['is_rep'], 'SampleName'].iloc[0]

    # clustermap
    g = sns.clustermap(df, mask=mask,
                       row_linkage=row_linkage,
                       col_linkage=col_linkage,

                       row_cluster=row_linkage is not None,

                       # yticklabels=False if len(sig_genes) == 0 else [g[:4] if (g[:6] != 'group_') & (g in sig_genes) else None for g in df.index],  # show only annotated gene labels
                       yticklabels=False if len(sig_genes) == 0 else [g[:4].replace('metE', '\nmetE').replace('cmpB', 'cmpB\n') if (g[:6] != 'group_') & (g[:3] not in ['rpl', 'rpm', 'rps']) & (g in sig_genes) else None for g in df.index],  # show only annotated gene labels
                       xticklabels=False if snps else ['*' if s == rep else None for s in df.columns],

                       row_colors=row_colors if species != 'Rep_449' else pd.DataFrame(['white' if g[:3] in ['rpl', 'rpm', 'rps'] else 'grey' for g in df.index], index=df.index, columns=['Ribosomal']),  # ribosomal proteins
                       col_colors=col_colors,

                       vmax=df.quantile(0.9).quantile(0.9) if snps & distance else None,
                       cmap=cmap,

                       cbar_kws={'label': 'Distance', 'orientation': 'horizontal', 'ticks': None} if distance else {})
    g.ax_heatmap.set_facecolor('black')
    g.ax_heatmap.tick_params(right=False, bottom=False)

    # dendrograms
    for ax, orientation in [(g.ax_col_dendrogram.axes, 'top'), (g.ax_row_dendrogram.axes, 'left')]:
        ax.clear()
        with plt.rc_context({'lines.linewidth': 0.5}):
            dendrogram(col_linkage if orientation == 'top' else row_linkage,
                       link_color_func=lambda x: link_colors[x] if distance | (orientation == 'top') else 'black',
                       ax=ax, orientation=orientation)
        ax.set_xticks([])
        ax.set_yticks([])

        if orientation == 'left':
            ax.invert_yaxis()

    # dimensions
    g.ax_heatmap.text(x=0.59, y=-0.035, s=f'n={df.shape[0]:,} {"" if distance else ("SNPs " if snps else "genes ")}x {df.shape[1]:,} {"samples" if snps else "strains"}', transform=g.ax_heatmap.transAxes)

    # title
    title = f'{segal_name(species)[0][3:]}'#{species}\n
    g.ax_col_dendrogram.set_title(title)

    # axes labels
    if snps:
        g.ax_heatmap.set_xlabel('Samples')
        g.ax_heatmap.set_ylabel('Samples' if distance else 'SNPs')
    else:
        g.ax_heatmap.set_xlabel('Strains', labelpad=-15)
        g.ax_heatmap.set_ylabel('Strains' if distance else 'Shell genes', labelpad=-35 if species != 'Rep_449' else 0)

    # position of plot
    plt.subplots_adjust(left=0.02, right=0.76, bottom=0.11, top=0.85)

    # position of distance legend
    if distance:
        g.cax.set_position([.26, .06, .5, .01])
    else:
        g.cax.set_visible(False)

    # significance
    if col_colors is not None:
        short_labels = []
        xlim = g.ax_col_colors.get_xlim()[1]
        max_len = max([len(text.get_text().split(' ~ ')[0]) for text in g.ax_col_colors.get_yticklabels()])
        for text in g.ax_col_colors.get_yticklabels():
            if '~' in text.get_text():
                label, tests = text.get_text().split(' ~ ')
                short_labels.append(label)
                colors_ = tests.split(' ')[::3]
                pvals = [float(p.split('=')[-1]) for p in tests.split(' ')[2::3]]
                y = text.get_position()[1]
                i = 0
                for c, p in zip(colors_, pvals):
                    if p < alpha:
                        # g.ax_col_colors.text(x=xlim, y=y+0.5, s=f'{" " * (5 + max_len + i*2)}*',
                        #                      color=c, fontsize=16)
                        i = i + 1
            else:
                label = text.get_text()
                short_labels.append(label)
        g.ax_col_colors.set_yticklabels(short_labels)

    if row_colors is not None:
        short_labels = []
        ylim = g.ax_row_colors.get_ylim()[0]
        max_len = max([len(text.get_text().split(' ~ ')[0]) for text in g.ax_row_colors.get_xticklabels()])
        for text in g.ax_row_colors.get_xticklabels():
            if '~' in text.get_text():
                label, tests = text.get_text().split(' ~ ')
                short_labels.append(label)
                colors_ = tests.split(' ')[::3]
                pvals = [float(p.split('=')[-1]) for p in tests.split(' ')[2::3]]
                x = text.get_position()[0]
                i = 0
                for c, p in zip(colors_, pvals):
                    if p < alpha:
                        # g.ax_row_colors.text(x=x-0.5, y=ylim, s=f'*{" " * (10 + max_len + i*2)}',
                        #                      color=c, fontsize=16, rotation='vertical', va='top')
                        i = i + 1
            else:
                label = text.get_text()
                short_labels.append(label)
        g.ax_row_colors.set_xticklabels(short_labels)

    g.ax_col_dendrogram.legend(handles=[
        Patch(facecolor='gold', edgecolor='black', label='Have gene'),
        Patch(facecolor='white', edgecolor='black', label='Missing gene'),
        Patch(facecolor='white', edgecolor='white', label=''),  # SPACE
        Patch(facecolor='royalblue', edgecolor='white', label='Old'),
        Patch(facecolor='lightsteelblue', edgecolor='white', label='Young'),
        Patch(facecolor='lightpink', edgecolor='white', label='Female'),
        Patch(facecolor='lightskyblue', edgecolor='white', label='Male'),
        Patch(facecolor='black', edgecolor='white', label='Unknown'),
    ], frameon=False, bbox_to_anchor=(-0.002, 1.1))

    subdir = 'snps' if snps else 'genes'
    subdir = (subdir + '_distance') if distance else subdir

    if not os.path.exists(os.path.join('figs', subdir)):
        os.makedirs(os.path.join('figs', subdir))
    plt.savefig(os.path.join('figs', subdir, species), pad_inches=0.5)

    if significant:
        if not os.path.exists(os.path.join('figs', subdir, 'significant')):
            os.makedirs(os.path.join('figs', subdir, 'significant'))
        plt.savefig(os.path.join('figs', subdir, 'significant', species), pad_inches=0.5)

    plt.close()

    return g.dendrogram_col.reordered_ind


def fig_snp_heatmap(
        # data
        species=None, snps=True, distance=True,
        # samples
        minimal_samples=100, maximal_imputation=0.1,
        # clustering
        metric='euclidean', method='average',
        # strain definition
        relative_change_threshold=0.05, min_samples_per_strain=50, max_samples_per_strain=0.8,
        # annotations
        col_annots=[], col_annots_covariates=[], col_annots_df=None, col_annots_cmaps={},
        row_annots=[], row_annots_covariates=[], row_annots_df=None, row_annots_cmaps={},
        ):

    col_dist, mask = get_distance(species, metric='diss' if snps else 'mash', maximal_imputation=maximal_imputation)
    if col_dist.shape[0] < minimal_samples:
        return {}
    col_linkage = get_linkage(col_dist, method)
    col_link_colors, color_names = find_strains(col_dist, col_linkage, relative_change_threshold, min_samples_per_strain, max_samples_per_strain)
    sample_colors = {i: (s, col_link_colors[np.where(col_linkage[:, :2] == i)[0][0] + col_dist.shape[0]]) for i, s in enumerate(col_dist.index)}

    if distance:
        data = col_dist
        row_linkage = col_linkage
    else:
        data = get_data(species, snps)
        data[list(set(col_dist.index) - set(data.columns))] = np.nan
        data = data[col_dist.index]
        data, mask = filter_data(data, sample_colors, min_samples_per_strain, '')
        if data.shape[0] == 0:
            return {}
        elif data.shape[0] == 1:
            row_linkage = None
        else:
            row_dist, mask = get_distance(species, metric, maximal_imputation, data)
            row_linkage = get_linkage(row_dist, method)

    pvalues = []
    expanded_row_annots = []
    for a in row_annots:
        pvalue, label = get_statistics(a, row_annots_covariates, row_annots_df, sample_colors, min_samples_per_strain)
        expanded_row_annots.append(label)
        # pvalues.append(pvalue)
    print(expanded_row_annots)
    if len(row_annots) > 0:
        row_colors = get_colors(row_annots_df.loc[row_annots_df.index.intersection(data.index), row_annots], row_annots_cmaps)
        row_colors.columns = expanded_row_annots
    else:
        row_colors = None

    sig_genes = []
    expanded_col_annots = []
    for a in col_annots:
        if a == 'Study':
            pvalue, label = 1, a
        else:
            pvalue, label = get_statistics(a, col_annots_covariates, col_annots_df, sample_colors, min_samples_per_strain)
        expanded_col_annots.append(label)
        pvalues.append(pvalue)
        if pvalue < alpha:
            new_sig_genes = filter_data(df=data,
                                    sample_colors={i: (v[0], v[1]) for i, v in col_annots_df[a].reset_index().replace('Unknown', np.nan).iterrows()},
                                    min_samples_per_strain=1,
                                    suffix=f'_{a}')
            if not ((species == 'Rep_449') & (a == 'Sex')):
                sig_genes = new_sig_genes[0].index.union(sig_genes)
        print(label)
    # print(expanded_col_annots)
    if len(col_annots) > 0:
        col_colors = get_colors(col_annots_df.loc[col_annots_df.index.intersection(data.columns), col_annots], col_annots_cmaps)
        col_colors.columns = expanded_col_annots
    else:
        col_colors = None

    reordered_ind = plot(df=data, mask=mask, distance=distance, snps=snps,
                         significant=(np.array(pvalues) < alpha).any(),  sig_genes=sig_genes,
                         row_linkage=row_linkage, col_linkage=col_linkage, link_colors=col_link_colors,
                         row_colors=row_colors, col_colors=col_colors)

    sample_colors = {sample_colors[i][0]: color_names[sample_colors[i][1]] for i in reordered_ind}

    return sample_colors


if __name__ == '__main__':
    snps = False
    distance = False
    study = None

    col_annots = row_annots = []
    col_annots_cmaps = row_annots_cmaps = {}
    col_annots_df = row_annots_df = None
    col_annots_covariates = row_annots_covariates = []

    # do once
    if snps:
        os.chdir(f'/home/saarsh/Analysis/strains/{study}')

        alpha = 0.05 / 102 if study == '10K' else 0.05 / 112

        with pd.HDFStore(os.path.join('data_frames', 'mb_dists.h5'), 'r') as hdf:
            all_species = [key[1:] for key in hdf.keys()]
        del hdf
        strains = pd.read_pickle(os.path.join('data_frames', 'meta.df'))[[]]

        pheno = pd.read_pickle(os.path.join('data_frames', 'pheno.df'))
        col_annots = col_annots_covariates = ['Age', 'Sex', 'BMI']
        col_annots_df = pheno[['age', 'sex', 'bmi']].rename(columns={'age': 'Age', 'sex': 'Sex', 'bmi': 'BMI'})

        if distance:
            row_annots = ['Has genome', 'Variable positions', 'Relative abundance'] if study == '10K' else ['Variable positions', 'Relative abundance']
            row_annots_cmaps = {'Variable positions': 'Greys', 'Relative abundance': 'Greys'}

            HG = pd.read_hdf(os.path.join(os.getcwd(), 'data_frames', 'has_genome.h5')) if study == '10K' else None
            SV = pd.read_hdf(strain_variability)
            RA = pd.read_hdf(os.path.join(os.getcwd(), 'data_frames', 'abundance.h5'))

        else:
            all_species = glob.glob(os.path.join('figs', 'snps_distance', 'significant', '*.png'))
            all_species = [s.split('/')[-1].split('.')[0] for s in all_species]

    else:
        os.chdir(f'/home/saarsh/Analysis/pangenomes')

        alpha = 0.05 / 191

        all_species = glob.glob(gene_presence_absence.format('*'))
        all_species = [file.split('/')[-3] for file in all_species]
        strains = pd.DataFrame(index=assemblies['SampleName'].unique())

        col_annots = ['Age', 'Sex', 'Study']
        assemblies['Sex'] = assemblies['Gender'].replace({'F': 'Female', 'M': 'Male'})
        assemblies['Study'] = assemblies['Source'].astype(str) + '_' + assemblies['StudyTypeID'].astype(str)
        studies = {s: i for i, s in enumerate(sorted(assemblies['Study'].unique()))}
        assemblies['Study'] = assemblies['Study'].replace(studies)

        # all_species = glob.glob(os.path.join('figs', 'genes_blue', 'significant', '*png'))
        # all_species = [file.split('/')[-1].split('.')[0] for file in all_species]
        all_species = ['Rep_449']#'Rep_721',

    # run
    for species in all_species:

        print(species)

        # do per species
        if snps:
            if distance:
                if study == '10K':
                    row_annots_df = HG[species].to_frame('Has genome').join(
                                    SV[species].to_frame('Variable positions'), how='outer').join(
                                    RA[species].to_frame('Relative abundance'), how='outer')
                else:
                    row_annots_df = SV[species].to_frame('Variable positions').join(
                                    RA[species].to_frame('Relative abundance'), how='outer')
        else:
            col_annots_df = \
            assemblies.loc[assemblies['cluster_s'] == int(species.split('_')[-1])].set_index('SampleName')[['Age', 'Sex', 'Study']]

        # plot
        sample_colors = fig_snp_heatmap(
            # data
            species=species, snps=snps, distance=distance,
            # samples
            minimal_samples=200, maximal_imputation=0.1,
            # clustering
            metric='euclidean', method='average',
            # strain definition
            relative_change_threshold=0.01, min_samples_per_strain=50, max_samples_per_strain=0.8,
            # annotations
            col_annots=col_annots, col_annots_covariates=col_annots_covariates, col_annots_df=col_annots_df,
            col_annots_cmaps={'Age': 'Blues', 'BMI': 'Reds', 'Study': 'nipy_spectral'},

            row_annots=row_annots, row_annots_covariates=row_annots_covariates, row_annots_df=row_annots_df,
            row_annots_cmaps=row_annots_cmaps,
        )

        # data frame
        if len(sample_colors) > 0:
            strains.loc[list(sample_colors.keys()), species] = list(sample_colors.values())

    if distance:
        strains.to_pickle(os.path.join('data_frames', 'strains.df'))
