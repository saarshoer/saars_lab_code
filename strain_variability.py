import os
import glob
import numpy as np
import pandas as pd

# PNP3
# counts_dir = '/net/mraid08/export/genie/LabData/Data/MBPipeline/PNP3_rerun_segata/MBSNP/CountsB/'
# species = 'SGB'

# Everything else
counts_dir = '/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/CountsS/'
species = 'Rep'

counts_fname = f'mb_snp_counts_{species}_*.h5'
summary_fname = 'mb_snb_strain_variability_R{}.h5'

min_reads = 3

summary_df = pd.DataFrame()
species_files = glob.glob(os.path.join(counts_dir, counts_fname))

for i, species_fname in enumerate(species_files):
    print(f'{i}/{len(species_files)}')
    variable_pos = total_pos = pd.Series(dtype=np.float64)
    with pd.HDFStore(species_fname, 'r') as hdf:

        for contig_part in hdf.keys():
            contig_df = hdf[contig_part]
            contig_df = contig_df[contig_df.sum(axis=1) >= min_reads]
            variable_pos = variable_pos.add(
                pd.DataFrame(contig_df.groupby('SampleName').apply(lambda s: ((s != 0).sum(axis=1) > 1).sum())),
                fill_value=0)
            total_pos = total_pos.add(
                pd.DataFrame(contig_df.groupby('SampleName').apply(lambda s: s.shape[0])),
                fill_value=0)

    summary_df = summary_df.join(variable_pos.div(total_pos).rename(
        columns={0: f'{species}_{species_fname.split("_")[-1].split(".")[0]}'}), how='outer')

    # if i == 1:
    #     break

summary_df.to_hdf(os.path.join(os.path.dirname(os.path.dirname(counts_dir)), summary_fname.format(min_reads)),
                  key='snps')

print('done')
