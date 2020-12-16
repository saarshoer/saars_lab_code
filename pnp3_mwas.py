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

    df_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/data_frames'.format(study_ids[0])

    def get_df(df_name, df_delta):
        df = pd.read_pickle(os.path.join(df_dir, df_name))
        df = df.xs(group, level='group') if group else df
        if df_delta:
            df = get_delta_df(df, '0months')
        else:
            df = df.xs(time_point, level='time_point') if time_point else df

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
    body = get_df('body.df', False)
    diet = get_df('diet.df', delta)
    diet.columns = [col.replace('%', '') for col in diet.columns]

    # mwas data frames
    joined_df = idx.join(blood).join(body).join(diet)
    joined_df = joined_df.reset_index().rename(columns={'person': 'RegistrationCode', 'time': 'Date'})\
        .set_index(['RegistrationCode', 'Date'])

    x = idx['sample'].tolist()
    y = prep_df(joined_df[y_cols])
    c = prep_df(joined_df[cov_cols])

    return x, y, c


class P:
    study_ids = ['PNP3']
    body_site = 'Gut'
    group = None
    time_point = '0months'
    delta = True

    y_cols = ['bt__hba1c']  # 'bmi' - participants were not suppose to change weight
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
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid


def gen_f(subjects_df, df):
    if subjects_df is not None:
        df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()
