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
    study_ids = ['PNP3']###not in new snps pipeline
    body_site = 'Gut'

    # dfs_dir = f'/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/strains/{study_ids[0]}/data_frames'

    # queue
    max_jobs = 250
    send_to_queue = True
    jobname = f'pnp3'
    verbose = False
    analyses_dir = config.analyses_dir
    work_dir_suffix = f'{"_".join(study_ids)}_diss'

    # species
    species_set = None#pd.read_pickle(os.path.join(dfs_dir, f'baseline_species.df')).index.tolist()
    # with pd.HDFStore('/net/mraid20/export/genie/LabData/Analyses/saarsh/20230129_113349/mb_dists.h5', 'r') as hdf:
    #     done_species = [s[1:] for s in hdf.keys()]
    # species_set = list(set(species_set) - set(done_species))
    # del hdf
    ignore_species = None
    filter_by_species_existence = True
    species_blocks = 1

    # subjects
    subjects_gen_f = lambda: SubjectLoader().get_data(study_ids=P.study_ids, groupby_reg='first')

    # samples
    # samples_set = pd.read_pickle(os.path.join(dfs_dir, 'meta.df')).index.tolist()
    samples_set = GutMBLoader().get_data('segal_species', study_ids=study_ids).df_metadata.index
    largest_sample_per_user = False
    min_positions_per_sample = 20000
    min_common_positions = 20000
    subsample_dir = ''#study_ids[0]
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

    # import glob
    # files = glob.glob('/net/mraid20/export/genie/LabData/Analyses/saarsh/20221221_131315/*.h5')
    # files = {i: f for i, f in enumerate(files)}
    # m._post_full_run(files)

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
