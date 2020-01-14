# needs to be run with the older version of python3
# py3_old /home/saarsh/Develop/Git/mine/mid_pipeline_snp.py
# maybe be removing the 3 from the following line in the LabMBPipeline.config_global
# python3 = globals().get('python3', '/usr/wisdom/python3/bin/python3')
# or something with the --use_general_python argument

import os
import pandas as pd

from LabUtils.Utils import mkdirifnotexists
from LabUtils.addloglevels import sethandlers
from LabMBPipeline import config_global as config
from LabMBPipeline.mmmbp import get_args, main
import logging

out_dir = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/IBD/mmmbp')
DFOut_dir = os.path.join(out_dir, 'DFOut_SNP')
tmp_dir = os.path.join(out_dir, 'tmp')

Email = 'saar.shoer@weizmann.com'

max_r = 550
General_params = ' --max_r {} --use_general_python '.format(max_r)

Modules = ' --module_seq "MID,UZP,SNP" '

if not os.path.exists(DFOut_dir):
    os.makedirs(DFOut_dir)

NUM_SAMPLES = 100

MID_params = '--mid_md_path ' + os.path.join(DFOut_dir, "PostUniteMetadata.df ") + \
             '--mid_input_path ' + os.path.join(out_dir, "tmp", "UNT", " ") + \
             '--mid_check_cont ' + \
             '--mid_ext .fastq.gz '

jobs_dir = os.path.join(tmp_dir, 'jobs')
mkdirifnotexists(jobs_dir)
os.chdir(jobs_dir)

md_df = pd.read_pickle(os.path.join(out_dir, 'DFOut', "PostUniteMetadata.df"))
md_df_done = pd.DataFrame(columns=md_df.columns)
md_df_done.to_pickle(os.path.join(DFOut_dir, "PostUniteMetadata.df"))

logf, warnstream = sethandlers(logging.INFO, True, True, file_prefix='mmmbp_')

with config.qp(jobname='mmmbp_', max_r=max_r, q=['himem7.q'],
               _python=config.python3, _mem_def='30G', _tryrerun=True) as q:
    q.startpermanentrun()

    cnt = 0
    while True:
        md_df = pd.read_pickle(os.path.join(out_dir, 'DFOut', "PostUniteMetadata.df"))
        md_df['RawReadLength'] = 100
        if os.path.exists(os.path.join(DFOut_dir, "PostUniteMetadata.df")):
            md_df_done = pd.read_pickle(os.path.join(DFOut_dir, "PostUniteMetadata.df"))
            md_df_done.to_pickle(os.path.join(DFOut_dir, 'PostUniteMetadata.df_{}'.format(cnt)))
            md_df_done['cnt'] = False

        inds = [i for i in md_df.index if i not in md_df_done.index][:NUM_SAMPLES]
        if len(inds) == 0:
            break
        cnt += 1
        print("{} old, {} new round {}".format(len(md_df_done), len(inds), cnt))

        md_df_done = pd.concat([md_df_done, md_df.loc[inds]], sort=False)
        md_df_done.loc[inds, 'cnt'] = True

        md_df_done.to_pickle(os.path.join(DFOut_dir, "PostUniteMetadata.df"))

        args, sm = \
            get_args((DFOut_dir + ' ' +
                      tmp_dir + ' ' +
                      Email +
                      General_params +
                      Modules +
                      MID_params
                      ).split())

        main(args, sm, q=q, job_name='mmmbp_mid_snp', logf=logf, warnstream=warnstream)

warnstream.close()
print("Done in %d rounds" % cnt)
