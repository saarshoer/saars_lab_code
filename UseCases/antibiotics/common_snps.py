import os
import glob
import numpy as np
import pandas as pd

run_type = 'within'
coef_col = 'Coef'

base_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics'

study1 = '10K'
study2 = 'Lifeline_deep'

dir1 = os.path.join(base_path, f'{study1}/{run_type}/raw_hdfs')
dir2 = os.path.join(base_path, f'{study2}/{run_type}/raw_hdfs')

files1 = [os.path.basename(f) for f in glob.glob(os.path.join(dir1, '*'))]
files2 = [os.path.basename(f) for f in glob.glob(os.path.join(dir2, '*'))]
common_files = set(files1) & set(files2)

dfs = []
for file in common_files:
    print(file)
    df1 = pd.read_hdf(os.path.join(dir1, file))[['N', coef_col]]
    df2 = pd.read_hdf(os.path.join(dir2, file))[['N', coef_col]]

    df = df1.join(df2, how='inner', lsuffix=f'_{study1}', rsuffix=f'_{study2}')
    dfs.append(df)

if run_type == 'within':
    pd.concat(dfs).to_hdf(os.path.join(base_path, f'common_snps_{run_type}.h5'), key='snps', complevel=9)
else:
    n_parts = 4
    for i in np.arange(n_parts):
        pd.concat(dfs[int(len(dfs)*i/n_parts):int(len(dfs)*(i+1)/n_parts)]).to_hdf(
            os.path.join(base_path, f'common_snps_{run_type}_pt{i+1}.h5'), key='snps', complevel=9)
