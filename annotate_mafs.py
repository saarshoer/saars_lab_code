import os
import glob
from LabQueue.qp import qp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import OralMBSNPLoader


def func():
    potential_species = glob.glob('/home/saarsh/Genie/LabData/Data/MBPipeline/Analyses/MBSNP/Oral/MAF/mb_snp_maf_SGB_*_R1_S100.h5')
    potential_species = ['SGB_' + s.split('_')[-3] for s in potential_species]
    done_species = glob.glob('/home/saarsh/Genie/LabData/Data/MBPipeline/Analyses/MBSNP/Oral/MAF/mb_snp_annot_maf_SGB_*_R1_S100.h5')
    done_species = ['SGB_' + s.split('_')[-3] for s in done_species]
    species = list(set(potential_species) - set(done_species))

    ld = OralMBSNPLoader()
    ld._gen_species_set_maf_annot_data(species, min_reads_per_snp=1, min_samples_per_snp_cached=100)
    # TODO: make sure the gene annotation loader is using the OralMBLoader and not the Gut


sethandlers(file_dir='/home/saarsh/Analysis/antibiotics/jobs/')
os.chdir('/home/saarsh/Analysis/antibiotics/jobs/')

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q']) as q:
    q.startpermanentrun()
    tkttores = {}
    tkttores[0] = q.method(func)
    for k, v in tkttores.items():
        q.waitforresult(v)
