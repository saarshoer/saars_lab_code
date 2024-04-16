import pandas as pd

species = 'SGB_6506'

vcf_path = f'/net/mraid20/export/genie/LabData/Data/MBPipeline/PNP3_rerun_segata/MBSNP/CountsB/mb_snp_counts_{species}.h5'
s1 = '946729_FD3159'
s2 = '800808_FD3093'

s1_counter = 0
shared_counter = 0

with pd.HDFStore(vcf_path, mode='r') as vcf:
    print(len(vcf.keys()))
    for key in vcf.keys():
        print(key)
        contig = vcf[key]
        contig = (contig[contig.sum(axis=1) > 2] > 0).any(axis=1).unstack('Position').fillna(False)
        if s1 in contig.index:
            s1_counter += contig.sum(axis=1).loc[s1]
            if s2 in contig.index:
                shared_counter += (contig.loc[s1] & contig.loc[s2]).sum()

print('s1', s1_counter)
print('shared', shared_counter)

diss_path = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/PNP3/data_frames/diss_20_corrected.df'
diss = pd.read_pickle(diss_path)

print(diss.xs(species, level='Species').xs(s1, level='SampleName1')['shared_pos'])

print('')
