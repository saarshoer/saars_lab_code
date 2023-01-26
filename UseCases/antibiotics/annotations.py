import os
import glob
import numpy as np
import pandas as pd
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataAnalyses.MBSNPs import mwas_annots

# parameters
base_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/10K/within'
maf_template = None#'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/Cache/global_mafs/10K/mb_snp_g_maf_{}_R1_S500.h5'

mwas_file_path = os.path.join(base_dir, 'mb_gwas_significant.h5')
# mwas_file_path = os.path.join(base_dir, 'raw_hdfs', 'mb_gwas_Rep_*_Rep_*.h5')

jobs_path = os.path.join(base_dir, 'jobs')


# run
os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='20G') as q:
    q.startpermanentrun()

    snps = q.method(mwas_annots.run, (mwas_file_path, base_dir, 'Gut', maf_template))
    q.waitforresult(snps)

# with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], _mem_def='50G') as q:#########if it fails again do one snp per geneid
#     jobs_path
#     q.startpermanentrun()
#
#     snps = {}
#     for file in glob.glob(mwas_file_path):
#         snps[file] = q.method(mwas_annots.run, (file, os.path.join(base_dir, 'tested_snps_annots_new'), 'Gut', maf_template))
#
#     for file in snps.keys():
#         q.waitforresult(snps[file])
#
# snps = glob.glob(os.path.join(base_dir, 'tested_snps_annots_new', 'mb_gwas_Rep_*_Rep_*.h5'))
# snps = [pd.read_hdf(file)[['GeneID', 'feature', 'gene', 'product', 'Preferred_name', 'best_og_cat']] for file in snps]
# snps = pd.concat(snps)
# ###snps = snps[[]].reset_index().drop('Y', axis=1).drop_duplicates(keep='first').assign(Y='Y').set_index(snps.index.names)
#
# snps['text'] = snps['gene'].fillna(snps['Preferred_name'].replace('-', np.nan)).fillna('    ').str[:4].str.replace('_', '').replace('    ', np.nan)
# snps['text'] = snps['text'].fillna(snps['product'].str.split(' ').str[0].replace('hypothetical', np.nan))
# snps = snps.dropna(subset=['text'])  # notice this drops un-annotated lines
#
# snps['best_og_cat'] = snps['best_og_cat'].replace('S', 'unknown').replace('-', 'unknown').fillna('unknown')
# snps.loc[(snps['best_og_cat'] != 'unknown') & (snps['best_og_cat'].str.len() > 1), 'best_og_cat'] = 'multiple'
#
# snps.to_hdf(os.path.join(base_dir, 'tested_snps_gene_annotations_very_short.h5'), key='snps', complevel=9)
