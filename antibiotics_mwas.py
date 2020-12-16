import os
import pandas as pd
from LabUtils.Utils import date2_dir
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS


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
    # max_samples_per_snp = 10000
    snp_set = None

    # covariates - required even if none
    covariate_gen_f = None
    # covariate_loaders = None
    # covariate_get_data_args = {}
    test_maf_cov_corr = False

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'countries': countries, 'groupby_reg': 'first'}

    output_cols = None
    # for cov in cov_cols:
    #     output_cols = output_cols + [cov + '_Pval', cov + '_Coef']


def y_gen_f_inner(subjects_df, y):
    # if subjects_df is not None:
    #     y = filter_dataframe(y, subjects_df)
    return LoaderData(pd.DataFrame(y), None)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()
