import os
import pandas as pd
from analysis import cat2binary
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
    study_ids = ['PNP3']
    countries = None

    # \\math17-lx\saarsh\Genie\LabData\Analyses\saarsh\PNP3_gut_mwas_6months_group
    df = pd.read_pickle(os.path.join(
        '/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/data_frames'.format(study_ids[0]),
        '{}_abundance.df'.format(body_site.lower())))
    df = df.xs('6months', level='time_point')
    x = list(df.index.get_level_values('sample').unique())
    y = df.reset_index().rename(columns={'person': 'RegistrationCode'}).set_index('RegistrationCode')['group']
    y = y.replace('mediterranean', 0).replace('algorithm', 1)

    # queue
    max_jobs = 300
    jobname = '{}_{}_mwas'.format('_'.join(study_ids), body_site.lower())
    send_to_queue = True#False
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    work_dir_suffix = jobname

    # species
    species_set = None#['SGB_1815', 'SGB_2318']
    ignore_species = None
    species_blocks = 5

    # samples
    samples_set = x
    largest_sample_per_user = False
    min_positions_per_sample = 0

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.98  # (Liron 0.99) Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 10  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 10  # (Liron 400)
    snp_set = None

    # covariates
    # covariate_loaders = None
    # covariate_get_data_args = {}
    test_maf_cov_corr = False  # necessary

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}

    # y
    y_gen_f = lambda subjects_df: y_gen_f_inner(subjects_df, P.y)
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid

    output_cols = None


def y_gen_f_inner(subjects_df, y):
    if subjects_df is not None:
        y = filter_dataframe(y, subjects_df)
    return LoaderData(pd.DataFrame(y), None)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()

    # M = MWASInterpreter(params=P, mwas_fname='mb_gwas.h5',
    #                     work_dir='/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_gut_mwas',
    #                     out_dir='/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/figs'.format(P.study_ids[0]),
    #                     mbsnp_loader=get_mbsnp_loader_class(P.body_site),
    #                     pval_col='Global_FDR', pval_cutoff=0.05,
    #                     SNPs_to_plot_dct={},
    #
    #                     do_manhattan_plot=True,
    #                     do_mafs_plot=False,
    #                     do_qq_plot=True,
    #                     do_volcano_plot=True,
    #
    #                     do_snp_annotations=False,
    #                     annotate_all_snps=False,
    #                     do_annotated_manhattan=False,
    #
    #                     get_extra_gene_info=False,
    #                     do_test_nonsynonymous_enrichment=False,
    #                     ).run()
