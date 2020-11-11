import os
# os.environ['OMP_NUM_THREADS'] = '2'  # Keep OMP from taking many CPUs.

import glob
import pandas as pd
from LabQueue.qp import qp
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.Loader import LoaderData
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.pandas_utils import filter_dataframe
from LabUtils.Utils import date2_dir, write_members
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class
from LabData.DataAnalyses.MBSNPs.MWASInterpreter import MWASInterpreter


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
    max_samples_per_snp = 10000
    snp_set = None

    # covariates
    # covariate_loaders = None
    # covariate_get_data_args = {}
    test_maf_cov_corr = False  # necessary

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'countries': countries, 'groupby_reg': 'first'}

    output_cols = None


def y_gen_f_inner(subjects_df, y):
    # if subjects_df is not None:
    #     y = filter_dataframe(y, subjects_df)
    return LoaderData(pd.DataFrame(y), None)


if __name__ == '__main__':
    # # creating the raw files
    # sethandlers(file_dir=config.log_dir)
    # m = MWAS(P)
    # work_dir = m.gen_mwas()

    # running the interpreter
    # TODO: perhaps change directories structure so to be one folder per species with within and between inside it
    run_type = 'between_species'

    input_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/anti_mwas_processed', run_type)
    output_path = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs', run_type)
    jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs/'

    def do(file_path):

        M = MWASInterpreter(params=P, mwas_fname=os.path.basename(file_path), work_dir=input_path, out_dir=output_path,
                            mbsnp_loader=get_mbsnp_loader_class(P.body_site),
                            pval_col='Global_Bonferroni', pval_cutoff=0.05,
                            SNPs_to_plot_dct={},

                            do_manhattan_plot=False,  # unnecessary given the annotated manhattan
                            # TODO: check if in between_species it results in seperate plots per x species
                            do_mafs_plot=False,  # broken
                            do_qq_plot=True,
                            do_volcano_plot=True,

                            do_snp_annotations=True,
                            annotate_all_snps=True,
                            do_annotated_manhattan=True,
                            # TODO: within plot does not look good - need to debug

                            get_extra_gene_info=True,
                            do_find_surrounding_genes=True,
                            do_test_nonsynonymous_enrichment=True,
                            ).run()

    os.chdir(jobs_path)
    sethandlers(file_dir=jobs_path)
    # TODO: between needs more memory
    with qp(jobname=run_type, _delete_csh_withnoerr=True, q=['himem7.q'], max_r=50) as q:
        q.startpermanentrun()
        tkttores = {}

        for file in glob.glob(os.path.join(input_path, 'SGB_9712.h5')):
            tkttores[file] = q.method(do, [file])

        for k, v in tkttores.items():
            q.waitforresult(v)

    print('done')
