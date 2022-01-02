import os
import pandas as pd
from analysis import get_delta_df
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPAnalyses

df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/data_frames'
df_suffix = '_20_corrected'
covariates = ['Age', 'Gender', 'Carbohydrates']


def gen_f(subjects_df, df_path):
    df = pd.read_pickle(df_path)
    # if subjects_df is not None:
    #     df = filter_dataframe(df, subjects_df)
    return LoaderData(df, None)


def is_y_valid(y, max_on_most_freq_val_in_col=0.99, min_on_non_freq_val_for_y=0.01, y_binary=None):#really effects the run time
    return True
    # major_count = y.value_counts().max()
    # minor_count = len(y) - major_count
    # return (major_count <= max_on_most_freq_val_in_col * len(y)) & (minor_count >= min_on_non_freq_val_for_y * len(y))


def get_data(body_site, time_point, delta, permute):

    def get_df(df_name, df_delta):
        df = pd.read_pickle(os.path.join(df_dir, df_name))

        if df_delta:
            df = get_delta_df(df, '0months', divide=True if df_delta == '/' else False)
        else:
            df = df.xs(time_point, level='time_point')

        names = ['person', 'Age', 'Gender'] if 'abundance' not in df_name else ['sample', 'person', 'Age', 'Gender']
        df = df.droplevel(df.index.names.difference(names))

        return df

    def prep_df(df):
        for col in df.columns:
            labels = sorted(df[col].dropna().unique().tolist())
            if len(labels) == 2:
                df[col] = df[col].replace(labels[0], 0).replace(labels[1], 1)
                print('column {} labels {} are now [0,1]'.format(col, labels))

        df.columns = [col.replace('% ', '') for col in df.columns]

        return df

    # raw data frames
    idx = get_df(f'{body_site.lower()}_species_abundance{df_suffix}.df', False)
    idx = idx.loc[:, (~idx.isna()).sum() >= 40]
    s = ['SGB_' + col.split('sSGB__')[-1] for col in idx.columns]
    idx = idx[[]]

    blood = get_df(f'blood{df_suffix}.df', delta)
    metabolites = get_df(f'metabolites{df_suffix}.df', delta)
    metabolites.columns = metabolites.columns.get_level_values('BIOCHEMICAL')
    cytokines = get_df(f'cytokines{df_suffix}.df', delta)
    cytokines.columns = cytokines.columns.get_level_values('ID')

    body1 = get_df(f'body{df_suffix}.df', False)  # for age and gender
    body2 = get_df(f'body{df_suffix}.df', delta)
    body = body1.reset_index(drop=False).set_index('person')[['Age', 'Gender']].join(body2)

    diet = get_df(f'diet{df_suffix}.df', delta)[['% Carbohydrates']]

    # mwas data frames
    joined_df = idx.join(blood).join(metabolites).join(cytokines).join(body).join(diet)
    joined_df = joined_df.reset_index(['sample']).rename(columns={'sample': 'SampleName'}).set_index('SampleName')
    joined_df.to_pickle(os.path.join(df_dir, 'mwas_perm', f'{body_site.lower()}_mwas_input.df'))
    joined_df = prep_df(joined_df)

    x = joined_df.index.tolist()
    if permute:
        org_idx = joined_df.index.copy()
        joined_df = joined_df.sample(frac=1, replace=False, random_state=1)
        joined_df.index = org_idx
    y = joined_df[joined_df.columns.difference(covariates)]
    c = joined_df[covariates]

    y.to_pickle(os.path.join(df_dir, 'mwas_perm', f'{body_site.lower()}_mwas_input_y.df'))
    c.to_pickle(os.path.join(df_dir, 'mwas_perm', f'{body_site.lower()}_mwas_input_c.df'))

    return x, y, c, s


class P:
    # data
    study_ids = ['PNP3']
    body_site = 'Gut'
    time_point = '0months'
    delta = True#'-'  # alternative is '/' or False
    permute = True

    output_cols = ['N', 'Coef', 'Pval', 'Coef_025', 'Coef_975']
    for cov in covariates:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']

    x, y, c, s = get_data(body_site, time_point, delta, permute)

    collect_data = False

    # queue
    max_jobs = 120
    send_to_queue = True
    analyses_dir = config.analyses_dir
    work_dir_suffix = f'{"_".join(study_ids)}_mwas_{body_site.lower()}'
    jobname = f'P{body_site}'
    verbose = False

    # fake world for special mwas
    body_site = study_ids[0]

    # species
    # done_species = pd.read_csv(
    #     '/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas_gut_0months/done_species.csv', index_col=0)
    # done_species = done_species.iloc[:, 0].to_list()
    species_set = s#list(set(s) - set(done_species))
    # print(f'running on {len(species_set)} new species')
    ignore_species = None
    filter_by_species_existence = False
    species_blocks = 1

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}

    # samples
    samples_set = x
    largest_sample_per_user = False
    min_positions_per_sample = 1  # None fails
    subsample_dir = ''
    other_samples_set = None  # ask eran

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 20  # in maf file name
    max_on_fraq_major_per_snp = 0.8  # Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 10  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 40  # (Liron 400)
    max_samples_per_snp = None
    snp_set = None

    # covariates
    covariate_gen_f = lambda subjects_df: gen_f(subjects_df, os.path.join(df_dir, 'mwas_perm', f'gut_mwas_input_c.df'))
    constant_covariate = False
    ret_cov_fields = True
    test_maf_cov_corr = False  # necessary

    # y
    y_gen_f = lambda subjects_df: gen_f(subjects_df, os.path.join(df_dir, 'mwas_perm', f'gut_mwas_input_y.df'))
    is_y_valid_f = is_y_valid  # Function that checks whether the analyzed y is valid
    max_on_most_freq_val_in_col = 0.8
    min_on_non_freq_val_for_y = 0.2

    # p-value
    max_pval_to_report = 0.05
    max_pval_to_detailed = None

    # others
    compute_pairwise_distances = False
    mwas_data_clusterer = None
    groupby_reg = None


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)

    # m = MWAS(P)
    # work_dir = m.gen_mwas()
    # print(work_dir)

    work_dir = f'/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas_gut_change_perm2'

    MBSNPAnalyses(P, work_dir).post_full_run_recovery_from_files()
    df = pd.read_hdf(os.path.join(work_dir, 'mb_gwas.h5'))
    df = df[df['Y_Bonferroni'] <= 0.05]
    df.to_hdf(os.path.join(work_dir, 'mb_gwas_significant.h5'), key='sig')

    # P.collect_data = True
    # P.snp_set = pd.read_hdf(os.path.join(work_dir, 'mb_gwas_significant.h5'))
    # P.work_dir_suffix = f'{P.work_dir_suffix}_data'
    # print(P.snp_set.index.get_level_values('Y').value_counts())
    # print(P.snp_set.index.get_level_values('Species').value_counts())
    # d = MWAS(P)
    # work_dir = d.gen_mwas()
    # print(work_dir)

    # sig_results = pd.read_hdf(os.path.join(work_dir, 'mb_gwas_significant.h5'))
    # P.y = P.y[sig_results.index.get_level_values('Y').unique().tolist()]
    # max_pval_to_report = 1
    # P.work_dir_suffix = f'{P.work_dir_suffix}_all_pvals_subset_y'
    # d = MWAS(P)
    # work_dir = d.gen_mwas()
    # print(work_dir)

    # work_dir = os.path.join(work_dir, 'hdfs_all_pvals_sig_y')
    # MBSNPAnalyses(P, work_dir).post_full_run_recovery_from_files()

    print('done')
