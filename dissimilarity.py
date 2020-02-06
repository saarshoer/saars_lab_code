import os
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabUtils.Utils import date2_dir, write_members, to_list
from LabData.DataMergers.MultiDataLoader import MultiDataLoader
from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances


class P:
    species_set = ['SGB_2304']#None
    send_to_queue = False
    largest_sample_per_user = False

    body_site = 'Gut'
    min_subjects_per_snp = 400
    min_reads_per_snp = 3

    # body_site = 'Oral'
    # min_subjects_per_snp = 10
    # min_reads_per_snp = 1

    subjects_loaders = ['SubjectLoader']

    study_ids = ['PNP3']
    subjects_get_data_args = {
        'study_ids': study_ids,
    }


def data_gen(loaders, subjects_df=None, **kwargs):
    accepts_subjects_df = all([l != 'SubjectLoader' for l in to_list(loaders)])
    return MultiDataLoader(loaders, subjects_df=subjects_df, **kwargs).get_data() if accepts_subjects_df \
        else MultiDataLoader(loaders, **kwargs).get_data()


def gen_pairwise_dists():
    work_dir = os.path.join(config.analyses_dir, date2_dir())
    write_members(os.path.join(work_dir, 'PARAMS.txt'), P)

    subjects_gen_f = lambda: data_gen(P.subjects_loaders, **P.subjects_get_data_args)

    MBSNPPairwiseDistances(work_dir=work_dir,
                           send_to_queue=P.send_to_queue,

                           body_site=P.body_site,
                           min_reads_per_snp=P.min_reads_per_snp,
                           min_subjects_per_snp=P.min_subjects_per_snp,
                           largest_sample_per_user=P.largest_sample_per_user). \
        run(subjects_gen_f=subjects_gen_f, species_set=P.species_set)


if __name__ == '__main__':
    sethandlers(file_dir=config.log_dir)
    gen_pairwise_dists()
