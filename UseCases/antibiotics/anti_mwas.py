import os
import glob
import pandas as pd
from LabUtils.Utils import date2_dir
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.pandas_utils import filter_dataframe

df = pd.read_pickle('/home/saarsh/Analysis/antibiotics/URA/dl.df')
df_metadata = pd.read_pickle('/home/saarsh/Analysis/antibiotics/URA/dl.df_metadata')


def gen_f(subjects_df, df):
    # if subjects_df is not None:
    #     df = filter_dataframe(df, subjects_df)
    return LoaderData(pd.DataFrame(df), None)


def is_y_valid(y, max_on_most_freq_val_in_col=0.95):
    return y.value_counts().max() <= max_on_most_freq_val_in_col * len(y)


class P:# run each time with one x and its corresponding abundance (if exist in URA!!! - if not what???) and perhaps already skip his y

    def __init__(self, species):
        self.species = species

        # general
        self.body_site = 'Gut'
        self.study_ids = ['D2']
        self.countries = ['IL']

        self.cov_cols = ['age', 'gender', 'abundance']

        self.y = df[['SGB_10068']]  # only within
        # self.y = df[[col != species for col in df.columns]]  # only between
        # self.c = df_metadata[cov_cols].dropna(how='any')
        self.c = df_metadata.join(df[[self.species]].rename(columns={self.species: 'abundance'}), how='inner')[self.cov_cols].dropna(how='any')
        # TODO: fix because it currently takes the pray instead of predator abundance
        # queue
        self.max_jobs = 500  # so to take no more than half the cluster's memory
        self.jobname = 'anti_mwas'
        self.send_to_queue = True#False
        self.work_dir = os.path.join(config.analyses_dir, 'anti_mwas_abundance', self.species)
        self.work_dir_suffix = self.jobname

        # species
        self.species_set = None#SGB_14399-1.61GB(smallest), SGB_4866-4.54GB, SGB_1815-50GB
        self.ignore_species = None
        self.filter_by_species_existence = False
        self.species_blocks = 1

        # subjects
        self.subjects_loaders = ['SubjectLoader']
        self.subjects_get_data_args = {'study_ids': self.study_ids, 'countries': self.countries, 'groupby_reg': 'first'}

        # samples
        self.samples_set = self.y.index.tolist()
        self.largest_sample_per_user = True
        self.min_positions_per_sample = 0

        # SNPs
        self.min_reads_per_snp = 1  # in maf file name
        self.min_subjects_per_snp_cached = 500  # in maf file name
        self.max_on_fraq_major_per_snp = 0.98  # (Eran 0.98, Liron 0.99) Max fraction of major allele frequency in analyzed samples
        self.min_on_minor_per_snp = 100  # (Liron 50) Min number of analyzed samples with a minor allele
        self.min_subjects_per_snp = 1000  # (Liron 400)
        # self.max_samples_per_snp = 10000
        self.snp_set = None

        # covariates - required even if none
        if self.cov_cols:
            self.covariate_gen_f = lambda subjects_df: gen_f(subjects_df, self.c)
        else:
            self.covariate_gen_f = None
        # self.covariate_loaders = None
        # self.covariate_get_data_args = {}
        self.test_maf_cov_corr = False  # necessary

        # y
        self.y_gen_f = lambda subjects_df: gen_f(subjects_df, self.y)
        self.is_y_valid_f = is_y_valid  # Function that checks whether the analyzed y is valid
        self.max_on_most_freq_val_in_col = 0.95  # make sure it has the same value as in is_y_valid_f

        # output
        self.output_cols = ['N', 'Coef', 'Pval', 'Coef_025', 'Coef_975']
        for cov in self.cov_cols:
            self.output_cols = self.output_cols + [cov + '_Pval', cov + '_Coef']


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)

    maf_species = glob.glob('/net/mraid08/export/genie/LabData/Data/MBPipeline/Analyses/MBSNP/Gut/MAF1/mb_snp_maf_SGB_*_R1_S500.h5')
    maf_species = [s.split('maf_')[1].split('_R1')[0] for s in maf_species]
    ura_species = df.columns
    for s in set(maf_species) & set(ura_species):
        # if s != 'SGB_10068':
        m = MWAS(P(species=s))
        work_dir = m.gen_mwas()
