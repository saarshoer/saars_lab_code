import os
import pandas as pd
from LabUtils.Utils import mkdirifnotexists
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.GutMBLoader import GutMBLoader
from LabData.DataLoaders.OralMBLoader import OralMBLoader

from LabMBPipeline import config_global as config
from LabMBPipeline.mmmbp import get_args, main

import logging

# touch
body_site = 'Gut'
run_name = 'Hanson'
MAX_JOBS = 250  # this is jobs, not threads


# do not touch
loader = GutMBLoader if body_site == 'Gut' else OralMBLoader

post = '{}/DFOut/PostUniteMetadata.df' if body_site == 'Gut' else 'Oral/{}/DFOut/PostUniteMetadata.df_no_URB'
post = os.path.join(config.mbpipeline_data_base, post.format(run_name))

read_len = {'Hanson': 75, 'NovaSeq': 100,  # Gut
            'Nextseq': 75, 'Novaseq': 100}  # Oral
read_len = read_len[run_name]

urb_num_mapped_to_subsample = {'Hanson': 8000000, 'NovaSeq': 5000000,  # Gut
                               'Nextseq': 4000000, 'Novaseq': 2500000}  # Oral
urb_num_mapped_to_subsample = urb_num_mapped_to_subsample[run_name]

urb_min_mapped_to_retain = {'Hanson': 1000000, 'NovaSeq': 1000000,  # Gut
                            'Nextseq': 500000, 'Novaseq': 500000}  # Oral
urb_min_mapped_to_retain = urb_min_mapped_to_retain[run_name]

base_dir = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/pipeline/', body_site, run_name)


# prep
mkdirifnotexists(base_dir)
mkdirifnotexists(os.path.join(base_dir, 'tmp2'))
mkdirifnotexists(os.path.join(base_dir, 'tmp2', 'jobs'))
mkdirifnotexists(os.path.join(base_dir, 'DFOut'))

df = pd.read_pickle(post)
samples = loader().get_data('segata_species', study_ids=['PNP3']).df.index
if body_site == 'Gut':
    oral_samples = OralMBLoader().get_data('segata_species', study_ids=['PNP3']).df.index
    samples = set(samples) - set(oral_samples)
df = df.loc[df.index.isin(samples), ['URB' not in col and 'SNPB' not in col and 'SNB' not in col for col in df.columns]]
if df['RawReadLength'].unique() != read_len:
    print('all hell broke loose')
    1/0
df['cnt'] = True
df.to_pickle(os.path.join(base_dir, 'DFOut', 'PostUniteMetadata.df'))
print(f'Running on {df.shape[0]}')


# pipeline
Email = ' saar.shoer@weizmann.ac.il'
General_params = f' --max_r {MAX_JOBS} --use_general_python '
Modules = ' --module_seq "MID,UZP,URB,SNB" '

os.chdir(os.path.join(base_dir, 'tmp2', 'jobs'))
logf, warnstream = sethandlers(logging.INFO, True, True, file_prefix='mmmbp_')

with config.qp(jobname=run_name, max_r=MAX_JOBS, q=['himem7.q'], _tryrerun=True, delay_batch=10) as q:
    q.startpermanentrun()

    MID_params = '--mid_md_path ' + os.path.join(base_dir, 'DFOut', 'PostUniteMetadata.df ') + \
                 '--mid_input_path ' + os.path.join(os.path.dirname(os.path.dirname(post)), 'tmp2', 'UNT', ' ') + \
                 '--mid_ext .fastq.gz ' + \
                 '--mid_check_cont '

    URB_params = f' --urb_num_mapped_to_subsample {urb_num_mapped_to_subsample} ' \
                 f' --urb_min_mapped_to_retain {urb_min_mapped_to_retain} ' \
                 ' --urb_run_type LargeOrNewGenusSGBs '

    SNB_params = ' --snb_run_type LargeOrNewGenusSGBs '

    args, sm = \
        get_args((os.path.join(base_dir, 'DFOut', ' ') +
                  os.path.join(base_dir, 'tmp2', ' ') +
                  Email +
                  General_params +
                  Modules +
                  MID_params +
                  URB_params +
                  SNB_params
                  ).split())

    print("max q jobs %d" % q._qp__max_r)
    main(args, sm, q=q, logf=logf, warnstream=warnstream)

warnstream.close()

print("Done")
print()
