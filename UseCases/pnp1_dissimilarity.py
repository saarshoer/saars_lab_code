import os
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabUtils.Utils import date2_dir, write_members, to_list

study_ids = ['PNP1']
species_set = None

subjects_loaders = ['SubjectLoader']
subjects_get_data_args = {'study_ids': study_ids, 'groupby_reg': 'first'}


class P:

    body_site = 'Gut'
    send_to_queue = True#False
    largest_sample_per_user = True

    min_reads_per_snp = 3

    min_common_positions = 20000
    min_positions_per_sample = 20000

    work_dir = os.path.join(config.analyses_dir, date2_dir())


def data_gen(loaders, subjects_df=None, **kwargs):
    from LabData.DataMergers.MultiDataLoader import MultiDataLoader
    accepts_subjects_df = all([l != 'SubjectLoader' for l in to_list(loaders)])
    return MultiDataLoader(loaders, subjects_df=subjects_df, **kwargs).get_data() if accepts_subjects_df \
        else MultiDataLoader(loaders, **kwargs).get_data()


def gen_pairwise_dists():
    write_members(os.path.join(P.work_dir, 'PARAMS.txt'), P)

    subjects_gen_f = lambda: data_gen(subjects_loaders, **subjects_get_data_args)

    from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances
    MBSNPPairwiseDistances(**dict((key, value) for key, value in P.__dict__.items()
                                  if not key.startswith('__')))\
        .run(subjects_gen_f=subjects_gen_f, species_set=species_set)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    gen_pairwise_dists()
