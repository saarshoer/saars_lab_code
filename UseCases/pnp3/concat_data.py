import os
import glob
import pandas as pd

base_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/PNP3_mwas_oral_small_data'
print(base_path)

hdfs = []
for file in glob.glob(os.path.join(base_path, 'mb_gwas_data_tmp*_N1.h5')):
    hdfs.append(pd.read_hdf(file))

pd.concat(hdfs).to_hdf(os.path.join(base_path, 'mb_gwas_data.h5'), key='snps')