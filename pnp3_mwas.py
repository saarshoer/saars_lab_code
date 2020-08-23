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
    body_site = 'Oral'
    study_ids = ['PNP3']
    countries = None

    time_point = '0months'
    label = 'group'

    # \\math17-lx\saarsh\Genie\LabData\Analyses\saarsh\PNP3_{oral/gut}_mwas_{0/6}months_group
    df = pd.read_pickle(os.path.join(
        '/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/data_frames'.format(study_ids[0]),
        '{}_abundance.df'.format(body_site.lower())))
    df = df.xs(time_point, level='time_point')
    x = list(df.index.get_level_values('sample').unique())
    y = df.reset_index().rename(columns={'person': 'RegistrationCode'}).set_index('RegistrationCode')[label]
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
    min_subjects_per_snp_cached = 500 if body_site == 'Gut' else 100  # in maf file name
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
    # sethandlers(file_dir=config.log_dir)
    # m = MWAS(P)
    # work_dir = m.gen_mwas()

    folder = '{}_{}_mwas_{}_{}'.format(P.study_ids[0], P.body_site.lower(), P.time_point, P.label)
    M = MWASInterpreter(params=P, mwas_fname='mb_gwas.h5',
                        work_dir=os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/', folder),
                        out_dir=os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/figs/'
                                             .format(P.study_ids[0]), folder),
                        mbsnp_loader=get_mbsnp_loader_class(P.body_site),
                        pval_col='Pval', pval_cutoff=0.05,
                        SNPs_to_plot_dct={},

                        do_manhattan_plot=True,
                        do_mafs_plot=False,  # broken
                        do_qq_plot=True,
                        do_volcano_plot=True,

                        do_snp_annotations=True,
                        annotate_all_snps=True,
                        do_annotated_manhattan=True,

                        get_extra_gene_info=True,
                        do_test_nonsynonymous_enrichment=True,
                        ).run()
