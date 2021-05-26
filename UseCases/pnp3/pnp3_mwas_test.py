import os
import numpy as np
import pandas as pd
from analysis import get_delta_df
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.pandas_utils import filter_dataframe
from LabUtils.timeutils import sample_date_19000101


def gen_f(subjects_df, df):
    if subjects_df is not None:
        df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


def is_y_valid(y, max_on_most_freq_val_in_col, min_on_non_freq_val_for_y, y_binary):
    return y.value_counts().max() <= max_on_most_freq_val_in_col * len(y)


def get_data(study_ids, body_site, group, time_point, y_cols, cov_cols, delta, log, sign, baseline):

    df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/data_frames'.format(study_ids[0])

    def get_df(df_name, df_delta):
        df = pd.read_pickle(os.path.join(df_dir, df_name))
        df = df.xs(group, level='group') if group else df
        baseline_df = df.xs('0months', level='time_point')

        if df_delta:
            if df_name != 'clinical_paper.df':
                divide = df_delta == '/'
                df = get_delta_df(df, '0months', divide)
            else:
                if df_delta == '-':
                    df = df.xs('6months', level='time_point')
                    df['time_point'] = f'6months-0months'
                elif df_delta == '/':
                    M0 = df.xs('0months', level='time_point')
                    MT = df.xs('6months', level='time_point')
                    MT = M0[MT.isna() == False].add(MT[M0.isna() == False])
                    df = MT[M0.isna() == False].div(M0[MT.isna() == False])
                    df['time_point'] = f'6months/0months'
                df = df.set_index('time_point', append=True)
            df = np.log2(df).replace(np.inf, np.nan).replace(-np.inf, np.nan) if log else df
            df = (df > 0).replace(True, 1).replace(False, 0).dropna() if sign else df
            df = 100*df.div(baseline_df) if baseline else df
        else:
            if df_name != 'clinical_paper.df' or (df_name == 'clinical_paper.df' and time_point == '0months'):
                df = df.xs(time_point, level='time_point') if time_point else df
            else:
                M0 = df.xs('0months', level='time_point')
                MT = df.xs(time_point, level='time_point')
                df = M0[MT.isna() == False].add(MT[M0.isna() == False])

        return df

    def prep_df(df):
        for col in df.columns:
            labels = df[col].unique()  # TODO: .dropna(), sort
            if len(labels) == 2:
                df[col] = df[col].replace(labels[0], 0).replace(labels[1], 1)
                print('column {} labels {} are now [0,1]'.format(col, labels))

        return df

    # raw data frames
    idx = get_df('{}_abundance.df'.format(body_site.lower()), False)
    indices2keep = ['person', 'group'] if delta else ['person', 'time_point', 'group']
    idx = idx[[]].reset_index(list(set(idx.index.names) - set(indices2keep)))
    blood = get_df('blood.df', delta)
    body1 = get_df('body.df', False)  # for age and gender
    body2 = get_df('body.df', delta)
    body = body1[['age', 'gender']].join(body2[list(set(body2.columns) - set(['age', 'gender']))])
    diet = get_df('diet.df', delta)
    clinical = get_df('clinical_paper.df', delta)  # data here is already in delta - special handling inside function

    # mwas data frames
    joined_df = idx.join(blood).join(body).join(diet).join(clinical)
    if body_site == 'Oral':
        joined_df['time'] = sample_date_19000101()
    joined_df = joined_df.reset_index().rename(columns={'person': 'RegistrationCode', 'time': 'Date'})\
        .set_index(['RegistrationCode', 'Date'])
    joined_df.columns = [col.replace('%', '') for col in joined_df.columns]

    x = idx['sample'].tolist()
    y = prep_df(joined_df[y_cols])
    c = prep_df(joined_df[cov_cols])

    return x, y, c


class P:
    # data
    study_ids = ['PNP3']
    body_site = 'Gut'
    group = None
    time_point = '0months'
    delta = '/'  # alternative is '/' or False
    log = False
    sign = False
    baseline = False

    # y_cols = \
    #     ['bt__neutrophils_', 'bt__hdl_cholesterol', 'bt__potassium',
    #      'bt__ldl_cholesterol', 'bt__lymphocytes_', 'bt__platelets',
    #      'bt__sodium', 'bt__hct', 'bt__mch', 'bt__mcv', 'bt__tsh',
    #      'bt__triglycerides', 'bt__eosinophils_', 'bt__rdw', 'bt__creatinine',
    #      'bt__ast_got', 'bt__mchc', 'bt__albumin', 'bt__crp_synthetic',
    #      'bt__glucose', 'bt__fructosamine', 'bt__monocytes_', 'bt__insulin',
    #      'bt__alt_gpt', 'bt__hemoglobin', 'bt__crp_hs', 'bt__total_cholesterol',
    #      'bt__rbc', 'bt__mean_platelet_volume', 'bt__hba1c', 'bt__wbc',
    #      'bt__basophils_', 'bt__calcium'] + \
    #     ['BMI', 'BMR', 'BP_dia', 'BP_sys', 'Fat_perc', 'HOMAIR', 'HbA1C_CGM', 'HeartRate',
    #      'Hips', 'OGTT', 'Time_above_140', 'Tot_Chol_to_HDL_Ratio', 'Waist', 'Weight']
    y_cols = ['Time_above_140']

    cov_cols = ['age', 'gender', 'carbohydrates']

    output_cols = ['N', 'Coef', 'Pval', 'Coef_025', 'Coef_975']
    for cov in cov_cols:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']

    x, y, c = get_data(study_ids, body_site, group, time_point, y_cols, cov_cols, delta, log, sign, baseline)

    collect_data = False
    max_pval_to_report = 1

    # queue
    max_jobs = 50
    send_to_queue = False
    work_dir = config.analyses_dir
    description = 'regular' if not delta else \
                 (('subtraction' if not sign else 'sign') if delta != '/' else \
                 ('division' if not log else 'log2division'))
    work_dir_suffix = f'{"_".join(study_ids)}_mwas_{body_site.lower()}_{time_point}_{description}'
    jobname = f'{body_site[0]}{description}'

    # species
    species_set = ['SGB_15342']
    ignore_species = None
    filter_by_species_existence = True
    species_blocks = 1

    # samples
    samples_set = x
    largest_sample_per_user = True
    min_positions_per_sample = 1
    subsample_dir = ''

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500 if body_site == 'Gut' else 100  # in maf file name
    max_on_fraq_major_per_snp = 0.99  # Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 10  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 10  # (Liron 400)
    max_samples_per_snp = None
    snp_set = pd.DataFrame(columns=['Y', 'Species', 'Contig', 'Position', 'fake col'])
    snp_set.loc[0] = ['Time_above_140', 'SGB_15342', 'C_3_P2', 135112, 'fake value']
    snp_set = snp_set.set_index(['Y', 'Species', 'Contig', 'Position'])

    # covariates
    covariate_gen_f = lambda subjects_df: gen_f(subjects_df, P.c)
    test_maf_cov_corr = False  # necessary

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}

    # y
    y_gen_f = lambda subjects_df: gen_f(subjects_df, P.y)
    is_y_valid_f = is_y_valid  # Function that checks whether the analyzed y is valid
    max_on_most_freq_val_in_col = 0.99

    # ?
    min_on_non_freq_val_for_y = None
    groupby_reg = None
    mwas_data_clusterer = None
    verbose = False
    constant_covariate = True
    max_pval_to_detailed = 1
    ret_cov_fields = True


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)

    m = MWAS(P)
    work_dir = m.gen_mwas()
    print(work_dir)

    # expected result
    # Y
    # Species
    # Contig
    # Position
    # N
    # Coef
    # Pval
    # Coef_025
    # Coef_975
    # age_Pval
    # age_Coef
    # gender_Pval
    # gender_Coef
    # carbohydrates_Pval
    # carbohydrates_Coef
    # Global_Bonferroni
    # Global_FDR
    # Time_above_140
    # SGB_15342
    # C_3_P2
    # 135112
    # 57 - 0.006522358
    # 8.55E-08 - 0.008625591 - 0.004419125
    # 0.544780841 - 0.005799249
    # 0.851511081
    # 0.028729071
    # 0.033428877
    # 0.65860491
    # 8.55E-08
    # 8.55E-08

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
