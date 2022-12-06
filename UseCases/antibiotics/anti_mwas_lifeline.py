import os
import pandas as pd
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from anti_mwas_functions import gen_cov_f, gen_y_f

study = 'Lifeline_deep'
df_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/data_frames'#/permuted'  # for both types of permutations


class P:

    snps = pd.read_pickle(os.path.join(df_dir, 'snps.df'))

    # general
    body_site = 'Gut'
    study_ids = [study]

    within = False#########don't forget to change the mem_def accordingly
    cov_cols = ['CONSTANT', 'age', 'gender'] if within else ['CONSTANT', 'age', 'gender', 'abundance']
    permute = 'permute' in df_dir

    countries = None
    collect_data = False

    # queue
    max_jobs = 250
    verbose = False
    send_to_queue = True
    jobname = 'anti_mwas_within' if within else 'anti_mwas_between'
    jobname = jobname + '_permuted' if permute else jobname
    jobname = study + '_' + jobname
    analyses_dir = config.analyses_dir
    work_dir_suffix = jobname

    # species
    species_set = snps.columns.tolist()
    ignore_species = None
    filter_by_species_existence = False
    species_blocks = 1  # most be 1 for this analysis, do not change!

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}

    # samples
    samples_set = snps.index.tolist()
    largest_sample_per_user = False
    min_positions_per_sample = None
    subsample_dir = study
    other_samples_set = None
    select_n_rand_samples = None

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.95  # Max fraction of major AND minor allele frequency in analyzed samples
    min_on_minor_per_snp = 50  # Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 500
    snp_set = pd.read_hdf(os.path.join(os.path.dirname(df_dir), 'within' if within else 'between', 'mb_gwas_significant.h5'))[[]] if collect_data else None
    # snp_set = snp_set.loc[snp_set.index.get_level_values('Species').isin([])]  # collect data make up
    # snp_set = snp_set.groupby('Species').apply(lambda data: data.iloc[:int(data.shape[0]/2)]).droplevel(0)
    # validation snps
    # snp_set = pd.read_hdf(os.path.join(os.path.dirname(os.path.dirname(df_dir)), '10K', 'within' if within else 'between', 'mb_gwas_significant_validation.h5'))
    # snp_set = snp_set.loc[snp_set['validation_level'] == 'tested', []]

    # covariates
    covariate_gen_f = lambda species: gen_cov_f(df_dir, species, P.within)
    constant_covariate = True
    ret_cov_fields = True
    test_maf_cov_corr = False  # necessary

    # y
    y_gen_f = lambda species: gen_y_f(df_dir, species, P.within)
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid
    max_on_most_freq_val_in_col = 0.95  # make sure it has the same value as in is_y_valid_f
    min_on_non_freq_val_for_y = 50

    # p-value
    max_pval_to_report = 10**-5 if permute else 1
    max_pval_to_detailed = None

    # others
    compute_pairwise_distances = False
    mwas_data_clusterer = None
    groupby_reg = None

    # output
    output_cols = ['N', 'R2', 'Pval', 'Coef', 'Coef_025', 'Coef_975']
    for cov in cov_cols:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']
        if 'abundance' in cov:
            output_cols = output_cols + [cov + '_Coef_025',  cov + '_Coef_975']

    del snps


if __name__ == '__main__':
    
    sethandlers(file_dir=config.log_dir)

    m = MWAS(P)
    m.gen_mwas()
