import os
import glob
from LabQueue.qp import qp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import OralMBSNPLoader


def func():
    species = glob.glob('/home/saarsh/Genie/LabData/Data/MBPipeline/Analyses/MBSNP/Oral/MAF/mb_snp_maf_SGB_*_R1_S100.h5')
    species = ['SGB_' + s.split('_')[-3] for s in species]

    ld = OralMBSNPLoader()
    ld._gen_species_set_maf_annot_data(species, min_reads_per_snp=1, min_samples_per_snp_cached=100)


sethandlers(file_dir='/home/saarsh/Analysis/antibiotics/jobs/')
os.chdir('/home/saarsh/Analysis/antibiotics/jobs/')

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q']) as q:
    q.startpermanentrun()
    tkttores = {}
    tkttores[0] = q.method(func)
    for k, v in tkttores.items():
        q.waitforresult(v)
