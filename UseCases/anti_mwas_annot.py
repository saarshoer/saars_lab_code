import os
import mwas_annot
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

# TODO: perhaps should be in jupyter or as part of analysis


def run():

    # parameters
    from UseCases.anti_mwas import P
    output_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed'

    x_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_raw/mb_gwas_SGB_143?.h5'
    y_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/*/SGB_*.h5'

    # run
    snps = mwas_annot.find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path)
    snps = mwas_annot.add_alleles(snps, P, output_path)  # early because it uses contig_with_part
    snps = mwas_annot.add_contig_without_part(snps, output_path)
    snps = mwas_annot.add_surrounding_genes(snps, output_path)
    snps = mwas_annot.flatten_current_genes(snps, output_path)
    snps = mwas_annot.add_annotations(snps, output_path)
    snps = mwas_annot.add_codons(snps, P, output_path)  # late because it depends on column 'feature' from annotations
    snps = mwas_annot.add_amino_acids(snps, output_path)

    return snps


jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs'

os.chdir(jobs_path)
sethandlers()

with fakeqp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='50G') as q:
    q.startpermanentrun()
    q.waitforresult(q.method(run))
