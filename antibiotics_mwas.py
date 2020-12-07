import os
# os.environ['OMP_NUM_THREADS'] = '2'  # Keep OMP from taking many CPUs.

import glob
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from matplotlib.lines import Line2D
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.pandas_utils import filter_dataframe
from LabUtils.Utils import date2_dir, write_members
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class
from LabData.DataAnalyses.MBSNPs.MWASInterpreter import MWASInterpreter


class P:
    # general
    body_site = 'Gut'
    study_ids = ['D2']
    countries = ['IL']

    # queue
    max_jobs = 100  # so to take no more than half the cluster's memory
    jobname = 'anti_mwas'
    send_to_queue = False#True
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    work_dir_suffix = jobname

    # species
    species_set = None#SGB_14399-1.61GB(smallest), SGB_4866-4.54GB, SGB_1815-50GB
    ignore_species = None
    # [file.split('mb_gwas_')[-1].split('.')[0]
    # for file in glob.glob('/home/saarsh/Genie/LabData/Analyses/saarsh/anti_mwas_raw/*h5')]

    species_blocks = 1

    # y
    y = pd.read_pickle('/home/saarsh/Analysis/antibiotics/URA/dl.df')
    y_gen_f = lambda subjects_df: y_gen_f_inner(subjects_df, P.y)
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid

    # samples
    samples_set = y.index.tolist()
    largest_sample_per_user = True
    min_positions_per_sample = 0

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.98  # (Liron 0.99) Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 100  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 1000  # (Liron 400)
    max_samples_per_snp = 10000
    snp_set = None

    # covariates
    # covariate_loaders = None
    # covariate_get_data_args = {}
    test_maf_cov_corr = False  # necessary

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'countries': countries, 'groupby_reg': 'first'}

    output_cols = None


def y_gen_f_inner(subjects_df, y):
    # if subjects_df is not None:
    #     y = filter_dataframe(y, subjects_df)
    return LoaderData(pd.DataFrame(y), None)


def text_func_annotated_between(df: pd.DataFrame):
    text = f"{df.shape[0]} SNPs in {df.index.get_level_values('Species').unique().shape[0]} species\n" + \
           f"N = between {df['N'].min()} and {df['N'].max()} samples per SNP\n" + \
           f"{sum(df['Global_Bonferroni'] <= 0.05)} SNPs passed 0.05 Global Bonferroni cutoff"
    return text


def text_func_annotated_within(df: pd.DataFrame):
    text = f"{sum(df['Global_Bonferroni'] <= 0.05)}/{df.shape[0]} SNPs passed 0.05 Global Bonferroni cutoff\n" + \
           f"N = between {df['N'].min()} and {df['N'].max()} samples per SNP"
    return text


def draw_func_annotated(sp_df: pd.DataFrame, **kwargs):
    d = dict()

    black = (0.0, 0.0, 0.0, 0.7)
    white = (1.0, 1.0, 1.0, 0.7)
    blue = (0.0, 0.0, 1.0, 0.7)
    red = (1.0, 0.0, 0.0, 0.7)

    marker_dict = {'intergenic': '*', 'synonymous': 'o', 'non-synonymous': 'o', 'non-protein': 'v'}
    facecolor_marker_dict = {'intergenic': black, 'synonymous': white, 'non-synonymous': black, 'non-protein': black}
    facecolor_dict = {True: red, False: blue}
    facecolor_label_dict = {True: 'positive relationship', False: 'negative relationship'}
    label_dict = {'intergenic': 'Intergenic (no annotated function)',
                  'synonymous': 'Protein coding - synonymous mutation',
                  'non-synonymous': 'Protein coding - non-synonymous mutation',
                  'non-protein': 'Non-protein coding (rRNA, tRNA etc.)'}

    d['marker'] = sp_df['annotation'].map(marker_dict).values
    d['s'] = -np.log10(sp_df[kwargs['pval_col']]).values
    # d['facecolor'] = sp_df['annotation'].map(facecolor_dict).values
    d['facecolor'] = (sp_df['Coef'] > 0).map(facecolor_dict).values
    d['edgecolor'] = black
    d['linewidths'] = 1
    d['alpha'] = 0.5

    legend_mutation = [Line2D([0], [0], linewidth=0, label=label_dict[k], alpha=d['alpha'],
                     markersize=10, marker=marker_dict[k], markeredgewidth=d['linewidths'],
                     markerfacecolor=facecolor_marker_dict[k], markeredgecolor=d['edgecolor']) for k in marker_dict.keys()]

    legend_coefficient = [Line2D([0], [0], linewidth=0, label=facecolor_label_dict[k], alpha=d['alpha'],
                     markersize=10, marker='s', markeredgewidth=0,
                     markerfacecolor=facecolor_dict[k], markeredgecolor=white) for k in facecolor_dict.keys()]

    legend = legend_mutation + legend_coefficient

    return sp_df, d, legend


if __name__ == '__main__':
    # # creating the raw files
    # sethandlers(file_dir=config.log_dir)
    # m = MWAS(P)
    # work_dir = m.gen_mwas()

    # running the interpreter
    run_type = 'between_species'

    input_path = os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed', run_type)
    output_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs', run_type)
    jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs/'

    # annotations
    annotations_df = pd.read_hdf('/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/snps_codon.h5')
    # annotation handling is done with contig without part while mwas data frame is contig with part
    annotations_df = annotations_df.rename(columns={'Contig_with_part': 'Contig'})
    annotations_df = annotations_df.reset_index('Contig', drop=True).set_index('Contig', append=True).\
        reorder_levels(['Species', 'Contig', 'Position'])
    # annotations that can be deduced (so not saved in the hdf files)
    annotations_df['NonSymMutation'] = annotations_df['MajorAA'] != annotations_df['MinorAA']
    annotations_df['annotation'] = annotations_df['NonSymMutation'].map({False: 'synonymous', True: 'non-synonymous'})
    annotations_df.loc[annotations_df['GeneDistance'] < 0, 'annotation'] = 'intergenic'
    annotations_df.loc[(annotations_df['GeneDistance'] >= 0) & (annotations_df['feature'] != 'CDS'), 'annotation'] = \
        'non-protein'  # (df['feature'] != 'CDS') does not cover all options
    # the only thing needed moving forward
    annotations_df = annotations_df['annotation']

    def do(file_path):

        M = MWASInterpreter(params=P, mbsnp_loader=get_mbsnp_loader_class(P.body_site),
                            mwas_fname=os.path.basename(file_path), work_dir=input_path,
                            out_dir=output_path, save_by_y=True,
                            pval_col='Global_Bonferroni', pval_cutoff=0.05,

                            do_manhattan_plot=True,
                            do_volcano_plot=True,
                            do_qq_plot=True,
                            # do_mafs_plot=False,  # broken

                            # do_snp_annotations=True,
                            # annotate_all_snps=True,
                            # do_annotated_manhattan=True,  # unnecessary given the modular manhattan

                            # get_extra_gene_info=True,
                            # do_find_surrounding_genes=True,
                            # do_test_nonsynonymous_enrichment=True,
                            )

        M._mwas_df = pd.read_hdf(os.path.join(M._work_dir, M._mwas_fname))
        M._mwas_df = M._mwas_df.join(annotations_df, on=['Species', 'Contig', 'Position'])
        M._draw_func = draw_func_annotated
        M._text_func = text_func_annotated_between if run_type == 'between_species' else text_func_annotated_within

        M.run()


    # queue

    os.chdir(jobs_path)
    sethandlers(file_dir=jobs_path)
    # TODO: between needs more memory
    with fakeqp(jobname=run_type, _delete_csh_withnoerr=True, q=['himem7.q'], max_r=50) as q:
        q.startpermanentrun()
        tkttores = {}

        for file in glob.glob(os.path.join(input_path, 'SGB_9712.h5')):#10068
            tkttores[file] = q.method(do, [file])

        for k, v in tkttores.items():
            q.waitforresult(v)

    print('done')
