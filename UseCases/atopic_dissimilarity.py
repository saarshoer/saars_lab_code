import os
import pandas as pd
from LabUtils.Utils import write_members
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.GutMBLoader import GutMBLoader
from LabData.DataLoaders.SubjectLoader import SubjectLoader
from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances


class P:
    # data
    # study_ids = ['AD_FMT2']  # first step
    # study_ids = ['AD_FMT']  # second step
    # study_ids = ['PNP1']  # control population
    # study_ids = ['AD_FMT', 'AD_FMT2']  # to compare donors or Rep_231
    body_site = 'Gut'

    # queue
    max_jobs = 250
    send_to_queue = True
    jobname = f'diss'
    verbose = False
    analyses_dir = config.analyses_dir
    work_dir_suffix = 'atopic'#f'{"_".join(study_ids)}_diss'

    # species
    species_set = None
    ignore_species = None
    filter_by_species_existence = True
    species_blocks = 1

    # subjects
    subjects_gen_f = lambda: SubjectLoader().get_data(groupby_reg='first')#study_ids

    # samples
    samples_set = GutMBLoader().get_data('segal_species', study_ids=['AD_FMT', 'AD_FMT2']).df.index.tolist() + \
                  pd.read_csv('/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/atopic/AD_FMT2/data_frames/metadata_10k.csv')['stool'].dropna().tolist()
    largest_sample_per_user = False
    min_positions_per_sample = 20000
    min_common_positions = 20000
    subsample_dir = ''
    other_samples_set = None
    select_n_rand_samples = None

    # SNPs
    min_reads_per_snp = 3
    min_subjects_per_snp = 1

    # irrelevant params
    min_subjects_per_snp_cached = None
    max_on_fraq_major_per_snp = None
    min_on_minor_per_snp = None
    snp_set = None

    # y
    is_y_valid_f = None
    max_on_most_freq_val_in_col = None
    min_on_non_freq_val_for_y = None

    # others
    mwas_data_clusterer = None
    groupby_reg = None


if __name__ == '__main__':

    sethandlers(file_dir=config.log_dir)

    m = MBSNPPairwiseDistances(P)
    write_members(os.path.join(m._work_dir, 'PARAMS.txt'), P)
    m.run(
        species_blocks=P.species_blocks,
        species_set=P.species_set,
        ignore_species=P.ignore_species,
        snp_set=P.snp_set,
        samples_set=P.samples_set,
        other_samples_set=P.other_samples_set,
        subjects_gen_f=P.subjects_gen_f,
        max_jobs=P.max_jobs,
        jobname=P.jobname)

    # any additional parameter needs to be specifically checked that it gets to its target function
