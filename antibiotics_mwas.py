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
    max_jobs = 850  # so to take no more than half the cluster's memory
    jobname = 'anti_mwas'
    send_to_queue = True#False
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    work_dir_suffix = jobname

    # species
    species_set = ['SGB_4866']#None#SGB_14399-1.61GB(smallest), SGB_4866-4.54GB, SGB_1815-50GB
    ignore_species = None
    species_blocks = 1

    # samples
    samples_set = None
    largest_sample_per_user = True
    min_positions_per_sample = 0

    # SNPs
    min_reads_per_snp = 1  # in maf file name
    min_subjects_per_snp_cached = 500  # in maf file name
    max_on_fraq_major_per_snp = 0.98  # (Liron 0.99) Max fraction of major allele frequency in analyzed samples
    min_on_minor_per_snp = 100  # (Liron 50) Min number of analyzed samples with a minor allele
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
                          'largest_sample_per_user': largest_sample_per_user, 'min_col_val': -4, 'min_col_present': None,
                          'min_col_present_frac': None, 'top_frac_present_col': None, 'take_log': True,
                          'convert_to_binary': False, 'groupby_reg': None}}
                          # eliminate in full run
                          # 'cols': ['k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_putredinis|fSGB__673|gSGB__1686|sSGB__2318',
                          #          'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_dorei|fSGB__595|gSGB__1353|sSGB__1815']}
    is_y_valid_f = None  # Function that checks whether the analyzed y is valid

    output_cols = None


if __name__ == '__main__':
    # sethandlers(file_dir=config.log_dir)
    # m = MWAS(P)
    # work_dir = m.gen_mwas()

    folder = '{}_{}_MAF'.format(P.jobname, '5GB' if P.species_set[0] is 'SGB_4866' else 'smallest')
    M = MWASInterpreter(params=P, mwas_fname='mb_gwas.h5',
                        work_dir=os.path.join('/net/mraid08/export/genie/LabData/Analyses/saarsh/', folder),
                        out_dir=os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/figs/'
                                             .format(P.study_ids[0]), folder),
                        mbsnp_loader=get_mbsnp_loader_class(P.body_site),
                        pval_col='Global_FDR', pval_cutoff=0.05,
                        SNPs_to_plot_dct={},

                        do_manhattan_plot=False,
                        do_mafs_plot=False,  # broken
                        do_qq_plot=False,
                        do_volcano_plot=True,

                        do_snp_annotations=False,
                        annotate_all_snps=False,
                        do_annotated_manhattan=True,

                        get_extra_gene_info=True,
                        do_test_nonsynonymous_enrichment=True,
                        ).run()
