import os
import pandas as pd
from analysis import get_delta_df
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS


df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/data_frames'
df_suffix = '_20_corrected'
covariates = ['Age', 'Gender', 'Carbohydrates']


def gen_f(subjects_df, df):
    # if subjects_df is not None:
    #     df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


def is_y_valid(y, max_on_most_freq_val_in_col=0.99, min_on_non_freq_val_for_y=0.01, y_binary=None):
    major_count = y.value_counts().max()
    minor_count = len(y) - major_count
    return (major_count <= max_on_most_freq_val_in_col * len(y)) & (minor_count >= min_on_non_freq_val_for_y * len(y))


def get_data(body_site, time_point, delta):

    def get_df(df_name, df_delta):
        df = pd.read_pickle(os.path.join(df_dir, df_name))

        if df_delta:
            df = get_delta_df(df, '0months', divide=True if df_delta == '/' else False)
            df = df.droplevel(df.index.names.difference(['person']))
        else:
            df = df.xs(time_point, level='time_point')

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
    joined_df.to_pickle(os.path.join(df_dir, f'{body_site.lower()}_mwas_input.df'))
    joined_df = prep_df(joined_df)

    x = joined_df.index.tolist()
    y = joined_df[joined_df.columns.difference(covariates)]
    c = joined_df[covariates]

    return x, y, c, s


class P:
    # data
    study_ids = ['PNP3']
    body_site = 'Oral'
    time_point = '0months'
    delta = '-'  # alternative is '/' or False

    output_cols = ['N', 'Coef', 'Pval', 'Coef_025', 'Coef_975']
    for cov in covariates:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']

    x, y, c, s = get_data(body_site, time_point, delta)

    collect_data = False

    # queue
    max_jobs = 200
    send_to_queue = True
    analyses_dir = config.analyses_dir
    work_dir_suffix = f'{"_".join(study_ids)}_mwas_{body_site.lower()}'
    jobname = f'{body_site}'
    verbose = False

    # fake world for special mwas
    body_site = study_ids[0]

    # species
    done_species = pd.read_csv(
        '/net/mraid08/export/genie/LabData/Analyses/saarsh/20210823_102934_PNP3_mwas_oral/finished.csv', index_col=0)
    done_species = done_species.iloc[:, 0].to_list()
    species_set = list(set(s) - set(done_species))
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
    max_on_fraq_major_per_snp = 0.95  # Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 10  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 10  # (Liron 400)
    max_samples_per_snp = None
    snp_set = None

    # covariates
    covariate_gen_f = lambda subjects_df: gen_f(subjects_df, P.c)
    constant_covariate = False
    ret_cov_fields = True
    test_maf_cov_corr = False  # necessary

    # y
    y_gen_f = lambda subjects_df: gen_f(subjects_df, P.y)
    is_y_valid_f = is_y_valid  # Function that checks whether the analyzed y is valid
    max_on_most_freq_val_in_col = 0.95
    min_on_non_freq_val_for_y = 0.05

    # p-value
    max_pval_to_report = 0.1
    max_pval_to_detailed = None

    # others
    compute_pairwise_distances = False
    mwas_data_clusterer = None
    groupby_reg = None


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)

    m = MWAS(P)
    work_dir = m.gen_mwas()
    print(work_dir)

    # work_dir = '/home/saarsh/Genie/LabData/Analyses/saarsh/???_PNP3_mwas_gut'
    # P.collect_data = True
    # P.snp_set = pd.read_hdf(os.path.join(work_dir, 'mb_gwas.h5'))
    # if not P.snp_set.empty:
    #     P.snp_set = P.snp_set[P.snp_set['Global_FDR'] <= 0.1]
    #     if not P.snp_set.empty:
    #         P.snp_set.to_hdf(os.path.join(work_dir, 'mb_gwas_significant.h5'), key='snps')
    #         print(P.snp_set.index.get_level_values('Y').value_counts())
    #         print(P.snp_set.index.get_level_values('Species').value_counts())
    #         P.work_dir_suffix = f'{P.work_dir_suffix}_data'
    #         d = MWAS(P)
    #         work_dir = d.gen_mwas()
