import os
import pandas as pd
from analysis import get_delta_df
from LabUtils.Utils import date2_dir
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.pandas_utils import filter_dataframe
from LabData.DataAnalyses.MBSNPs.MBSNPSpeciesSeparation import get_species_covariate, get_contig_covariate


def gen_f(subjects_df, df):
    if subjects_df is not None:
        df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


def is_y_valid(y, max_on_most_freq_val_in_col=0.99):
    return y.value_counts().max() <= max_on_most_freq_val_in_col * len(y)


def get_data(study_ids, body_site, group, time_point, y_cols, cov_cols, delta):

    # blood
    # ['bt__neutrophils_%', 'bt__hdl_cholesterol', 'bt__potassium',
    #  'bt__ldl_cholesterol', 'bt__lymphocytes_%', 'bt__platelets',
    #  'bt__sodium', 'bt__hct', 'bt__mch', 'bt__mcv', 'bt__tsh',
    #  'bt__triglycerides', 'bt__eosinophils_%', 'bt__rdw', 'bt__creatinine',
    #  'bt__ast_got', 'bt__mchc', 'bt__albumin', 'bt__crp_synthetic',
    #  'bt__glucose', 'bt__fructosamine', 'bt__monocytes_%', 'bt__insulin',
    #  'bt__alt_gpt', 'bt__hemoglobin', 'bt__crp_hs', 'bt__total_cholesterol',
    #  'bt__rbc', 'bt__mean_platelet_volume', 'bt__hba1c', 'bt__wbc',
    #  'bt__basophils_%', 'bt__calcium']

    # body
    # ['age', 'gender', 'sitting_blood_pressure_systolic', 'hips', 'body_fat',
    #  'bmi', 'bmr', 'height', 'sitting_blood_pressure_diastolic', 'waist',
    #  'trunk_fat', 'meetingtypeid', 'sitting_blood_pressure_pulse_rate',
    #  'weight']

    # diet
    # ['%proteins', '%lipids', '%carbohydrates']

    # clinical
    # ['ALT_GPT', 'AST_GOT', 'BMI', 'BMR', 'BP_dia', 'BP_sys',
    #  'Cholesterol___HDL', 'Cholesterol_total___blood', 'FPG', 'Fat_perc',
    #  'Fructosamine', 'HOMAIR', 'HbA1C_Blood', 'HbA1C_CGM', 'HeartRate',
    #  'Hips', 'Insulin', 'LDL_Cholesterol', 'OGTT', 'Time_above_140',
    #  'Tot_Chol_to_HDL_Ratio', 'Triglycerides', 'US', 'Waist', 'Weight']

    df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/data_frames'.format(study_ids[0])

    def get_df(df_name, df_delta):
        df = pd.read_pickle(os.path.join(df_dir, df_name))
        df = df.xs(group, level='group') if group else df
        if df_delta:
            df = get_delta_df(df, '0months')
        else:
            tp = time_point if df_name != 'clinical_paper.df' else '6months'
            df = df.xs(tp, level='time_point') if tp else df

        return df

    def prep_df(df):
        for col in df.columns:
            labels = df[col].unique()
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
    clinical = get_df('clinical_paper.df', False)  # data here is already in delta
    if delta:
        clinical['time_point'] = '6months-0months'
        clinical = clinical.set_index('time_point', append=True)

    # mwas data frames
    joined_df = idx.join(blood).join(body).join(diet).join(clinical)
    joined_df = joined_df.reset_index().rename(columns={'person': 'RegistrationCode', 'time': 'Date'})\
        .set_index(['RegistrationCode', 'Date'])
    joined_df.columns = [col.replace('%', '') for col in joined_df.columns]

    x = idx['sample'].tolist()
    y = prep_df(joined_df[y_cols])
    c = prep_df(joined_df[cov_cols])

    return x, y, c


class P:
    study_ids = ['PNP3']
    body_site = 'Gut'
    group = None
    time_point = '0months'  # does not apply to clinical
    delta = True  # for carbohydrates

    y_cols = \
        ['bt__neutrophils_', 'bt__hdl_cholesterol', 'bt__potassium',
         'bt__ldl_cholesterol', 'bt__lymphocytes_', 'bt__platelets',
         'bt__sodium', 'bt__hct', 'bt__mch', 'bt__mcv', 'bt__tsh',
         'bt__triglycerides', 'bt__eosinophils_', 'bt__rdw', 'bt__creatinine',
         'bt__ast_got', 'bt__mchc', 'bt__albumin', 'bt__crp_synthetic',
         'bt__glucose', 'bt__fructosamine', 'bt__monocytes_', 'bt__insulin',
         'bt__alt_gpt', 'bt__hemoglobin', 'bt__crp_hs', 'bt__total_cholesterol',
         'bt__rbc', 'bt__mean_platelet_volume', 'bt__hba1c', 'bt__wbc',
         'bt__basophils_', 'bt__calcium'] +\
        ['sitting_blood_pressure_systolic', 'hips', 'body_fat',
         'bmi', 'bmr', 'sitting_blood_pressure_diastolic', 'waist',
         'trunk_fat', 'sitting_blood_pressure_pulse_rate', 'weight'] + \
        ['proteins', 'lipids', 'carbohydrates']  # carbohydrates is intentional to see if the system works

    # Clinical
    # y_cols = \
    #     ['ALT_GPT', 'AST_GOT', 'BMI', 'BMR', 'BP_dia', 'BP_sys',
    #      'Cholesterol___HDL', 'Cholesterol_total___blood', 'FPG', 'Fat_perc',
    #      'Fructosamine', 'HOMAIR', 'HbA1C_Blood', 'HbA1C_CGM', 'HeartRate',
    #      'Hips', 'Insulin', 'LDL_Cholesterol', 'OGTT', 'Time_above_140',
    #      'Tot_Chol_to_HDL_Ratio', 'Triglycerides', 'US', 'Waist', 'Weight']  # carbohydrates is intentional to see if the system works


    cov_cols = ['age', 'gender', 'carbohydrates']

    # coverage
    species_specific_cov_f = get_species_covariate
    contig_specific_cov_f = get_contig_covariate

    output_cols = ['N', 'Coef', 'Pval', 'Coef_025', 'Coef_975',
                   'av_sample_score_Pval', 'av_sample_score_Coef',
                   'contig_av_sample_score_Pval', 'contig_av_sample_score_Coef']
    for cov in cov_cols:
        output_cols = output_cols + [cov + '_Pval', cov + '_Coef']

    x, y, c = get_data(study_ids, body_site, group, time_point, y_cols, cov_cols, delta)

    # queue
    max_jobs = 100
    jobname = '{}_{}_mwas'.format('_'.join(study_ids), body_site.lower())
    send_to_queue = True
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    work_dir_suffix = jobname

    # species
    species_set = None
    ignore_species = None
    filter_by_species_existence = True
    species_blocks = 1

    # samples
    samples_set = x
    largest_sample_per_user = False
    min_positions_per_sample = 0

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500 if body_site == 'Gut' else 100  # in maf file name
    max_on_fraq_major_per_snp = 0.99  # Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 10  # (Liron 50) Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 10  # (Liron 400)
    max_samples_per_snp = None
    snp_set = None

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


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()
