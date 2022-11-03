import os
import glob
import pandas as pd

run_type = 'between'
coef_col = 'Coef'

base_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas'

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
    df1 = pd.read_hdf(os.path.join(dir1, file))[[coef_col]]
    df2 = pd.read_hdf(os.path.join(dir2, file))[[coef_col]]

    df = df1.join(df2, how='inner', lsuffix=f'_{study1}', rsuffix=f'_{study2}')
    dfs.append(df)

if run_type == 'within':
    pd.concat(dfs).to_hdf(os.path.join(base_path, f'common_snps_{run_type}.h5'), key='snps', complevel=9)
else:
    pd.concat(dfs[:int(len(dfs)/2)]).to_hdf(os.path.join(base_path, f'common_snps_{run_type}_pt1.h5'), key='snps', complevel=9)
    pd.concat(dfs[int(len(dfs)/2):]).to_hdf(os.path.join(base_path, f'common_snps_{run_type}_pt2.h5'), key='snps', complevel=9)
