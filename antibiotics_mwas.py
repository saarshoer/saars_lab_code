import os
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs.MWAS import MWAS
from LabUtils.Utils import date2_dir, write_members
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class
from LabData.DataAnalyses.MBSNPs.MWASInterpreter import MWASInterpreter


class P:
    # general
    body_site = 'Gut'
    study_ids = ['D2']
    countries = ['IL']

    # queue
    max_jobs = 300
    jobname = 'anti'
    send_to_queue = True#False
    work_dir = os.path.join(config.analyses_dir, date2_dir())

    # species
    species_set = None#['SGB_2318', 'SGB_1815']
    ignore_species = None
    species_blocks = 5

    # samples
    largest_sample_per_user = True
    min_positions_per_sample = 0

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.98  # (Liron 0.99) Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 100  # (Liron 50)Min number of analyzed samples with a minor allele
    min_subjects_per_snp = 1000  # (Liron 400)
    snp_set = None

    # covariates
    # covariate_loaders = None
    # covariate_get_data_args = {}
    test_maf_cov_corr = False  # necessary

    # subjects
    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {'study_ids': study_ids, 'countries': countries, 'groupby_reg': 'first'}

    # y
    y_loaders = [body_site + 'MBLoader']
    y_get_data_args = {**{'df': 'segata_species', 'study_ids': study_ids, 'genotek_vals': None, 'min_hgf_fraction': None,
                       'min_reads': None, 'nextera': None, 'reg_ids': None, #'subjects_df': None,
                       'largest_sample_per_user': largest_sample_per_user, 'min_col_val': None, 'min_col_present': None,
                       'min_col_present_frac': None, 'top_frac_present_col': None, 'take_log': False,
                       'convert_to_binary': False, 'groupby_reg': None}}
                       # eliminate in full run
                       # 'cols': ['k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_putredinis|fSGB__673|gSGB__1686|sSGB__2318',
                       #          'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_dorei|fSGB__595|gSGB__1353|sSGB__1815']}
    #min_reads
    take_log = True
    min_col_val = -4  # in log space
    # min_col_present/_frac = 0.05 # Select columns that are present in at least min_col_present/_frac of the samples
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid

    output_cols = None


if __name__ == '__main__':
    write_members(os.path.join(P.work_dir, 'PARAMS.txt'), P)
    sethandlers(file_dir=config.log_dir)
    m = MWAS(P)
    work_dir = m.gen_mwas()

    # M = MWASInterpreter(params=P, mwas_fname='mb_gwas.h5',
    #                     work_dir='/net/mraid08/export/genie/LabData/Analyses/saarsh/antibiotics_mwas_test',
    #                     out_dir='/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics',
    #                     mbsnp_loader=get_mbsnp_loader_class(P.body_site),
    #                     pval_col='Global_FDR', pval_cutoff=0.05,
    #                     SNPs_to_plot_dct={},
    #
    #                     do_manhattan_plot=True,
    #                     do_mafs_plot=False,
    #                     do_qq_plot=True,
    #                     do_volcano_plot=True,
    #
    #                     do_snp_annotations=False,
    #                     annotate_all_snps=False,
    #                     do_annotated_manhattan=False,
    #
    #                     get_extra_gene_info=False,
    #                     do_test_nonsynonymous_enrichment=False,
    #                     ).run()
