import os
import pandas as pd
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS

df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/data_frames/permuted'

def gen_cov_f(species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'meta.df'))[['age', 'gender']]
    if not within:
        df2 = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))[species]
        df = df.\
            join(df2.rename(columns={species[0]: 'abundance'})).\
            join(df2.rename(columns={species[0]: 'MAF_abundance'}))
    return LoaderData(df, None)


def gen_y_f(species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))
    if within:
        df = df[species]
    else:
        if species[0] in df.columns:
            df = df.drop(species[0], axis=1)
    return LoaderData(df, None)


def is_y_valid(y, max_on_most_freq_val_in_col, min_on_non_freq_val_for_y):#really affects the run time
    count_most = y.value_counts().max()
    return (count_most <= max_on_most_freq_val_in_col * len(y)) & \
           (count_most >= (1 - max_on_most_freq_val_in_col) * len(y)) & \
           (len(y) - count_most >= min_on_non_freq_val_for_y)


class P:

    snps = pd.read_pickle(os.path.join(df_dir, 'snps.df'))

    # general
    body_site = 'Gut'
    study_ids = ['10K']

    within = False#########don't forget to change the mem_def accordingly
    cov_cols = ['age', 'gender'] if within else ['age', 'gender', 'abundance', 'MAF_abundance']

    countries = None
    collect_data = False

    # queue
    max_jobs = 300
    verbose = False
    send_to_queue = True
    jobname = 'perm_anti_mwas'
    analyses_dir = config.analyses_dir
    work_dir_suffix = jobname

    # species
    species_set = snps.columns.tolist()
    ignore_species = None
    filter_by_species_existence = False
    species_blocks = 1

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}

    # samples
    samples_set = snps.index.tolist()
    largest_sample_per_user = False
    min_positions_per_sample = None
    subsample_dir = '10K'
    other_samples_set = None
    select_n_rand_samples = None

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.95#also minor  # Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 50  # Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 500
    max_samples_per_snp = None
    snp_set = None

    # covariates
    covariate_gen_f = lambda species: gen_cov_f(species, P.within)
    constant_covariate = True
    ret_cov_fields = True
    test_maf_cov_corr = False  # necessary

    # y
    y_gen_f = lambda species: gen_y_f(species, P.within)
    is_y_valid_f = is_y_valid  # Function that checks whether the analyzed y is valid
    max_on_most_freq_val_in_col = 0.95  # make sure it has the same value as in is_y_valid_f
    min_on_non_freq_val_for_y = 50

    # p-value
    max_pval_to_report = 10**-5
    max_pval_to_detailed = None

    # others
    compute_pairwise_distances = False
    mwas_data_clusterer = None
    groupby_reg = None

    # output
    output_cols = ['N', 'Pval', 'Coef', 'Coef_025', 'Coef_975']
    for cov in cov_cols:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']
        if cov not in ['age', 'gender']:
            output_cols = output_cols + [cov + '_Coef_025',  cov + '_Coef_975']

    del snps


if __name__ == '__main__':
    
    sethandlers(file_dir=config.log_dir)

    m = MWAS(P)
    m.gen_mwas()
