import glob
from LabData.DataLoaders.MBSNPLoader import OralMBSNPLoader

species = glob.glob('/home/saarsh/Genie/LabData/Data/MBPipeline/Analyses/MBSNP/Oral/MAF/mb_snp_maf_SGB_*_R1_S100.h5')
species = ['SGB_' + s.split('_')[-3] for s in species]

ld = OralMBSNPLoader()
ld._gen_species_set_maf_annot_data(species, min_reads_per_snp=1, min_samples_per_snp_cached=100)
