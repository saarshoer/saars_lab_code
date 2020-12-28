import os
import glob
import mwas_plots
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabUtils.addloglevels import sethandlers

# TODO: perhaps should be in jupyter or as part of analysis

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
    run_type = 'between_species'

    input_path = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed', run_type)
    output_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs', run_type)
    jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs/'
    annotations_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/snps_maf_codon_annotations.h5'

    # annotations
    print('starting annotations')
    annotations_df = mutation_type(annotations_path)
    print('finished annotations')

    # queue
    os.chdir(jobs_path)
    sethandlers(file_dir=jobs_path)

    with qp(jobname=run_type, _delete_csh_withnoerr=True, q=['himem7.q'], _mem_def='10G', _tryrerun=True, _num_reruns=3) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for file in glob.glob(os.path.join(input_path, 'SGB_*.h5')):  # 9702, 10068
            kwargs = {'mwas_fname': file, 'annotations_df': annotations_df, 'out_dir': output_path,
                      'manhattan_draw_func': draw_func_annotated,
                      'manhattan_text_func': text_func_annotated_between if run_type == 'between_species' else text_func_annotated_within,
                      'pval_col': 'Pval', 'pval_cutoff': 0.05/26068850133}
            tkttores[file] = q.method(mwas_plots.run, kwargs=kwargs)
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
