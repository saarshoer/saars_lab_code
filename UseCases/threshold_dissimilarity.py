import os
import pandas as pd
from LabUtils.Utils import write_members
from LabData import config_global as config
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.SubjectLoader import SubjectLoader
from LabData.DataAnalyses.MBSNPs.MBSNPAnalyses import MBSNPPairwiseDistances

samples_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/dissimilarities/threshold_samples.df'


class P:
    # data
    study_ids = ['PNP1']# don't forget this is included in the atopic run VCFS/before_potential_donors/
    body_site = 'Gut'

    # queue
    max_jobs = 250
    send_to_queue = True
    jobname = f'diss'
    verbose = False
    analyses_dir = config.analyses_dir
    work_dir_suffix = f'{"_".join(study_ids)}_diss_{body_site.lower()}'

    # species
    species_set = None
    ignore_species = None
    filter_by_species_existence = False
    species_blocks = 1

    # subjects
    subjects_gen_f = lambda: SubjectLoader().get_data(study_ids=P.study_ids, groupby_reg='first')

    # samples
    samples_set = pd.read_pickle(samples_path).index.tolist()
    largest_sample_per_user = False
    min_positions_per_sample = 1000
    min_common_positions = 1000
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


# if __name__ == '__main__':
#
#     import matplotlib.pyplot as plt
#     from LabData.DataLoaders.GutMBLoader import GutMBLoader
#
#     samples = GutMBLoader().get_data('segal_species').df_metadata
#     samples = samples[samples['Date'].astype(str).str[:4].isin(['2014', '2015']) & (samples['RunType'] == 'NextSeq') &
#                       # next two lines are just to filter out odd samples
#                       (samples['URSMapUsed'] == 7 * 10 ** 6) &
#                       ~((samples['IsGenotek'] == 0) & (samples['PE'] == False) & (samples['Nextera'] == 'False'))]
#     samples = samples.groupby(['RegistrationCode', 'Date']).filter(lambda g: (g.shape[0] == 2) &
#                                                                              (set(g['IsGenotek']) == {0, 1}) &
#                                                                              (set(g['PE']) == {False, True}))
#     samples = samples.sort_values(['RegistrationCode', 'Date'])
#     samples.to_pickle(samples_path)
#
#     print(samples[['RegistrationCode', 'Date']])
#     print(samples[['IsGenotek', 'PE', 'Nextera']].value_counts())
#     bad_cols = []
#     for col in samples.columns:
#         l = len(set(samples[col]))
#         if (l > 1) & (l < 80):
#             try:
#                 plt.figure()
#                 plt.hist(samples[col])
#                 plt.title(col)
#             except:
#                 bad_cols.append(col)
#     print(bad_cols)

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
