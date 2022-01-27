import os
import pandas as pd
from LabUtils.Utils import write_members
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.SubjectLoader import SubjectLoader
from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances


class P:
    # data
    # study_ids = ['AD_FMT2']  # first step
    # study_ids = ['AD_FMT']  # second step
    # study_ids = ['PNP1']  # control population
    study_ids = ['AD_FMT', 'AD_FMT2']  # to compare donors or Rep_231
    body_site = 'Gut'

    # queue
    max_jobs = 250
    send_to_queue = True
    jobname = f'diss'
    verbose = False
    analyses_dir = config.analyses_dir
    work_dir_suffix = f'{"_".join(study_ids)}_diss_{body_site.lower()}'

    # species
    # species_set = pd.read_csv('~/reps_by_time.csv').iloc[:, 0]
    # species_set = [f'Rep_{s}' for s in species_set]
    species_set = None
    # species_set = ['Rep_231']
    ignore_species = None
    filter_by_species_existence = False
    species_blocks = 1

    # subjects
    subjects_gen_f = lambda: SubjectLoader().get_data(study_ids=P.study_ids, groupby_reg='first')

    # samples
    samples_set = None
    # samples_set = ['FMTAD040_v0_fullrun', 'FMTAD041_v0_fullrun', 'FMTAD042_v0_fullrun',
    #                'FMTAD043_v0_fullrun', 'FMTAD044_v0_fullrun', 'FMTAD045_v0_fullrun',
    #                'FMTAD046_v0_fullrun', 'FMTAD047_v0_fullrun', 'FMTAD048_v0_fullrun',
    #                'FMTAD049_v0_fullrun', 'FMTAD050_v0_fullrun', 'FMTAD051_v0_fullrun',
    #                'FMTAD052_v0_fullrun', 'FMTAD053_v0_fullrun', 'FMTAD054_v0_fullrun',
    #                'FMTAD055_v0_fullrun', 'fmtad121_v2_fullrun', 'fmtad122_v2_fullrun',
    #                'fmtad123_v2_fullrun', 'fmtad124_v2_fullrun', 'fmtad125_v2_fullrun',
    #                'fmtad126_v2_fullrun', 'fmtad127_v2_fullrun', 'fmtad128_v2_fullrun',
    #                'fmtad129_v2_fullrun', 'fmtad130_v2_fullrun', 'fmtad131_v2_fullrun',
    #                'fmtad132_v2_fullrun']  # to compare donors
    largest_sample_per_user = False
    min_positions_per_sample = 20000
    min_common_positions = 20000
    subsample_dir = ''
    other_samples_set = None

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
