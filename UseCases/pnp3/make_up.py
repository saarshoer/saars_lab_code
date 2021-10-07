import os
import glob
import pandas as pd

base_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/20210916_135718_PNP3_mwas_oral'
print(base_path)

done_species = []
for file in glob.glob(os.path.join(base_path, 'mb_gwas_tmp*_N1_counts.h5')):
    new_species = pd.read_hdf(file).index.get_level_values('Species').unique().tolist()
    print(new_species)
    done_species = done_species + new_species

pd.DataFrame(done_species).to_csv(os.path.join(base_path, 'done_species.csv'))