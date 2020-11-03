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


def get_data(study_ids, body_site, time_point, idx_cols, blood_cols, body_cols, diet_cols, cov_cols):

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

    # data data frames
    abundance = pd.read_pickle(os.path.join(df_dir, '{}_abundance.df'.format(body_site.lower())))
    abundance = abundance.xs(time_point, level='time_point')  # if time_point is not None else abundance
    abundance.index = abundance.index.set_names('RegistrationCode', level='person')

    blood = pd.read_pickle(os.path.join(df_dir, 'blood.df'))
    blood = blood.xs(time_point, level='time_point')  # if time_point is not None else blood
    blood.index = blood.index.set_names('RegistrationCode', level='person')

    body = pd.read_pickle(os.path.join(df_dir, 'body.df'))
    body = body.xs(time_point, level='time_point')  # if time_point is not None else body
    body.index = body.index.set_names('RegistrationCode', level='person')

    diet = pd.read_pickle(os.path.join(df_dir, 'diet.df'))
    diet = diet.xs(time_point, level='time_point')  # if time_point is not None else diet
    diet.index = diet.index.set_names('RegistrationCode', level='person')

    # mwas data frames
    idx_ = abundance[[]].reset_index(list(set(abundance.index.names)-set(['RegistrationCode'])))[idx_cols]
    blood_ = blood[blood_cols]
    body_ = body[body_cols]
    diet_ = diet[diet_cols]

    cov_ = body[cov_cols]

    x = list(abundance.index.get_level_values('sample').unique())

    y = idx_.join(blood_, how='outer').join(body_, how='outer').join(diet_, how='outer')
    y.index = y.index.droplevel(list(set(y.index.names) - set(['RegistrationCode'])))
    for col in y.columns:
        labels = y[col].unique()
        if len(labels) == 2:
            y[col] = y[col].replace(labels[0], 0).replace(labels[1], 1)
            print('column {} labels {} are now [0,1]'.format(col, labels))

    c = cov_
    cov_.index = cov_.index.droplevel(list(set(cov_.index.names) - set(['RegistrationCode'])))

    return x, y, c


class P:
    # general
    body_site = 'Oral'#
    time_point = '0months'#

    study_ids = ['PNP3']

    idx_columns = []
    blood_cols = ['bt__hba1c']#maybe separately
    body_cols = ['bmi']#maybe separately
    diet_cols = []

    cov_cols = ['age', 'gender']

    x, y, c = get_data(study_ids, body_site, time_point,
                       idx_columns, blood_cols, body_cols, diet_cols,
                       cov_cols)

    # queue
    max_jobs = 300
    jobname = '{}_{}_mwas'.format('_'.join(study_ids), body_site.lower())
    send_to_queue = False#True
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

    output_cols = None


def gen_f(subjects_df, df):
    if subjects_df is not None:
        df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


if __name__ == '__main__':

    # def _build_covariate_gen_f(self):
        # return self._params.covariate_gen_f if not hasattr(self._params, 'covariate_gen_f') else \
        #     lambda subjects_df: \
        #         self._data_gen(self._params.covariate_loaders, subjects_df=subjects_df, **self._params.covariate_get_data_args)

    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()

    # folder = '{}_{}_mwas_{}_{}'.format(P.study_ids[0], P.body_site.lower(), P.time_point, P.label)
    # M = MWASInterpreter(params=P, mwas_fname='mb_gwas.h5',
    #                     work_dir=os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/', folder),
    #                     out_dir=os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/{}/figs/'
    #                                          .format(P.study_ids[0]), folder),
    #                     mbsnp_loader=get_mbsnp_loader_class(P.body_site),
    #                     pval_col='Pval', pval_cutoff=0.05,
    #                     SNPs_to_plot_dct={},
    #
    #                     do_manhattan_plot=True,
    #                     do_mafs_plot=False,  # broken
    #                     do_qq_plot=True,
    #                     do_volcano_plot=True,
    #
    #                     do_snp_annotations=True,
    #                     annotate_all_snps=True,
    #                     do_annotated_manhattan=True,
    #
    #                     get_extra_gene_info=True,
    #                     do_test_nonsynonymous_enrichment=True,
    #                     ).run()
