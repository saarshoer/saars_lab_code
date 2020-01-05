import os
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabUtils.Utils import date2_dir, write_members, to_list


class p():

    min_positions_per_sample = 20 if config.DEBUG else 20000
    min_common_positions = 20000
    send_to_queue = False if config.DEBUG else True

    species_list = None

    species_set = ['SGB_1814'] if config.DEBUG else species_list

    study_ids = ['PNP1']

    subjects_loaders = ['SubjectLoader']
    subjects_get_data_args = {
        'study_ids': study_ids,
        'groupby_reg': 'first',
    }


def data_gen(loaders, subjects_df=None, **kwargs):
    from LabData.DataMergers.MultiDataLoader import MultiDataLoader
    accepts_subjects_df = all([l != 'SubjectLoader' for l in to_list(loaders)])
    return MultiDataLoader(loaders, subjects_df=subjects_df, **kwargs).get_data() if accepts_subjects_df \
        else MultiDataLoader(loaders, **kwargs).get_data()


def gen_pairwise_dists():
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    write_members(os.path.join(work_dir, 'PARAMS.txt'), p)

    subjects_gen_f = lambda: data_gen(p.subjects_loaders, **p.subjects_get_data_args)

    from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances

    MBSNPPairwiseDistances(min_common_positions=p.min_common_positions, send_to_queue=p.send_to_queue,
                           min_positions_per_sample=p.min_positions_per_sample,
                           work_dir=work_dir,
                           largest_sample_per_user=False).\
        run(subjects_gen_f=subjects_gen_f, species_set=p.species_set)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    gen_pairwise_dists()
