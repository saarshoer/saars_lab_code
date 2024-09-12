import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.GutMBLoader import GutMBLoader


###segal representative genomes
###pangenomes


# study_ids = ['NICU']
# output_dir = '/net/mraid20/export/jasmine/card/NICU'
# jobs_dir = os.path.join(output_dir, 'jobs')
# lookout_paths = ['/net/mraid20/export/genie/LabData/Data/MBPipeline/NebNext_new/tmp2/UNT/*sample*fastq*',
#                  '/net/mraid20/export/genie/LabData/Data/MBPipeline/NovaSeq/tmp2/UNT/*sample*fastq*']

# study_ids = ['10K']
# output_dir = '/net/mraid20/export/jasmine/card/10K'

# only has 200/600 samples, not worth the run
# study_ids = ['10K_young']
# output_dir = '/net/mraid20/export/jasmine/card/10K_young'

# study_ids = ['AD_FMT', 'AD_FMT2']
# output_dir = '/net/mraid20/export/jasmine/card/AD_FMTS'

# lookout_paths = ['/net/mraid20/export/genie/LabData/Data/MBPipeline/PipelineRerun/tmp2/HGF/*sample*fastq*',
#                  '/net/mraid20/export/genie/LabData/Data/MBPipeline/NebNext_updated/tmp2/UNT/*sample*fastq*']

# study_ids = ['Lifeline_deep', 'Lifeline_Stool']
# output_dir = '/net/mraid20/export/jasmine/card/Lifelines'
#
# lookout_paths = ['/net/mraid20/export/genie/LabData/Data/MBPipeline/Lifeline_75bp/tmp2/UNT/*sample*fastq*',
#                  '/net/mraid20/export/genie/LabData/Data/MBPipeline/Lifeline/tmp2/UNT/*sample*fastq*']

threads = 16
rgi_exe = '/home/saarsh/Develop/Git/rgi/rgi'

# done once on August 18th, 2024
# rgi load \
#   --card_json ~/Analysis/card/card/card.json \
#   --card_annotation ~/Analysis/card/card_database_v3.2.9.fasta \
#   --card_annotation_all_models ~/Analysis/card/card_database_v3.2.9_all.fasta \
#   --wildcard_annotation ~/Analysis/card/wildcard_database_v4.0.2.fasta \
#   --wildcard_annotation_all_models ~/Analysis/card/wildcard_database_v4.0.2_all.fasta \
#   --wildcard_index ~/Analysis/card/wildcard/index-for-model-sequences.txt \
#   --wildcard_version 4.0.2 \
#   --amr_kmers ~/Analysis/card/wildcard/all_amr_61mers.txt \
#   --kmer_database ~/Analysis/card/wildcard/61_kmer_db.json \
#   --kmer_size 61


def do(f, s):
    os.chdir(output_dir)
    os.system(f'python {rgi_exe} bwt -1 {f} -n {threads} -o {s} --clean --include_wildcard > /dev/null 2>&1')
    os.remove(f'{s}.allele_mapping_data.json')
    bams = glob.glob(f'{s}.sorted.length_*.bam*')
    for bam in bams:
        os.remove(bam)


if __name__ == '__main__':

    # queue
    jobs_dir = os.path.join(output_dir, 'jobs')
    os.makedirs(jobs_dir, exist_ok=True)
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    samples = GutMBLoader().get_data('segal_species', study_ids=study_ids).df_metadata.index.tolist()

    with qp(jobname='rgi', _mem_def='20G', _trds_def=threads, _tryrerun=False) as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        i = 0
        for sample in samples:
            files = []
            for lp in lookout_paths:
                files = files + glob.glob(lp.replace('sample', sample))
            if len(files) == 0:
                print(f'missing sample {sample}')
            elif len(files) > 1:
                print(f'duplicated sample {sample}')
            else:
                if not os.path.exists(os.path.join(output_dir, f'{sample}.gene_mapping_data.txt')):
                    removes = glob.glob(os.path.join(output_dir, f'{sample}*'))
                    for remove in removes:
                        os.remove(remove)
                    # print(sample)
                    tkttores[sample] = q.method(do, (files[0], sample))
                i = i + 1
        print('finished sending jobs')

        print(f'{i}/{len(samples)} ({i / len(samples) * 100}) found')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v, _assert_on_errors=False)  # see https://github.com/arpcard/rgi/issues/256
        print('finished waiting for jobs')

    print('done')


# in new pipeline number of samples

# 1001: Beilinson_Cardio 427
# 1005: Colorectal_cancer_recovered 12
# 1006: BreastCancerRecovered 97
# 1007: BRCA 115
# 1008: Endometriosis 109
# 1009: 10K_young 181
# 1010: 10K_T1D 14

# 15: Garvan 147
# 17: BreastCancer 224
# 20: IBD_Ichilov 550
# 38: Cancer_Wolfson 96
# 39: Parkinson 146
# 44: BreastCancerProfiling 180
# 45: PancreaticCancer 273


# not in new pipeline

# 1002: OvarianCancer
# 1003: Ichilov-crc-recovered
# 1004: Ichilov-BreastCancer-recovered

# 1000: D2

# 33: Lifeline_ibd not in loader
# 34: Lifeline_celiac not in loader

# 1: PNP1
# 2: PNP2
# 3: PNP3
# 4: ElinavAS
# 5: MS
# 6: T1D
# 7: GDM
# 8: Cardio
# 9: MDG
# 11: FattyLiver
# 14: NYU
# 16: Longitudinal
# 19: HealthyChildren
# 21: LVAD
# 22: CHF
# 25: IBDMDB
# 26: CFS_UK
# 27: CFS_UK_control
# 28: ICI_French
# 31: TwinsUK
# 35: COVID19
# 36: CRAFT
# 37: AUS_PID
# 40: Bariatric
# 41: Propionic
# 46: Glioblastoma
# 49: T2D
# 50: Bread
# 51: Cardio-Lymphoma
# 52: Cardio-heart transplantation
# 98: PRIMM
# 99: Press
# 101: LongReads
