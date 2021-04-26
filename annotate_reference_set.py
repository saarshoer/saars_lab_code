import os
import glob
import datetime
import pandas as pd
from Bio import SeqIO
from LabQueue.qp import qp, fakeqp
from LabUtils.Utils import shell_command
from LabUtils.addloglevels import sethandlers
from LabData.config_global import analyses_dir
from LabData.CentralizedGenerators.prokka import prokka

# in shell - create a conda environment
# conda create -n prokka_env -c conda-forge -c bioconda -c defaults prokka

# in shell - add the prokka and eggnog environments to your .cshrc PATH, if not in file itself it will not work for jobs
# setenv PATH ~/Develop/anaconda_33/envs/prokk_env/bin:/net/mraid08/export/genie/Bin/eggnogv2/eggnog-mapper:${PATH}

# in shell - run this script with your regular python (the one that has LabQueue and LabUtils configured)
# /usr/wisdom/python3/bin/python ~/Develop/Git/LabData/LabData/CentralizedGenerators/annotate_reference_genome.py &

# in python -
reference_fasta = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/mash/Mar21_final/BuildURA/Rep_all/single_fastas'  # Segal
# reference_fasta = '/net/mraid08/export/jafar/Microbiome/Data/Databases/URA_IndicesAndScores/LargeOrNewGenusSGBs/SegataIndexLargeOrNewGenusSGBs.fa'  # Segata
# reference_fasta = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/URA_DBs/SegataAll/other_single_fasta'  # Segata all
output_dir = analyses_dir

jobs_dir = os.path.join(output_dir, 'jobs')

prokka_env = '/home/saarsh/Develop/anaconda_33/envs/prokk_env/bin'
eggnog_env = '/net/mraid08/export/genie/Bin/eggnogv2/eggnog-mapper'
prokka_cmd = 'prokka --cpus 2 --quiet --force --outdir {} --prefix prokka {} --rfam'
eggnog_cmd = 'emapper.py -i {} --output {} --dbmem --go_evidence all -m diamond --cpu 16'
prokka_version = 'prokka_1.14.6'
eggnog_version = 'eggnog_2.0.4'

# create fasta file per species
if os.path.isfile(reference_fasta):
    prev_rep = ''
    reference_file = open(reference_fasta, 'r')
    for species_object in SeqIO.parse(reference_file, 'fasta'):
        rep = species_object.id.split('_c')[0]
        new_rep = rep != prev_rep
        species_dir = os.path.join(output_dir, rep)
        if new_rep:
            os.mkdir(species_dir)
        species_path = os.path.join(species_dir, rep) + '.fa'
        mode = 'w' if new_rep else 'a'
        species_file = open(species_path, mode)
        species_file.write(f'>{"_".join(species_object.id.split("_", 4)[:4])}\n{species_object.seq}\n')
        species_file.close()
        prev_rep = rep
    reference_file.close()
elif os.path.isdir(reference_fasta):
    for rep in glob.glob(os.path.join(reference_fasta, '*')):
        species_dir = os.path.join(output_dir, os.path.basename(rep).replace('.fa', ''))
        os.mkdir(species_dir)
        os.chdir(species_dir)
        shell_command(f"ln -s {rep} {os.path.basename(rep)}")
        # in case you need to edit contig names
        # reference_file = open(rep, 'r')
        # species_file = open(os.path.basename(rep), 'w')
        # for i, species_object in enumerate(SeqIO.parse(reference_file, 'fasta')):
        #     species_file.write(f">{os.path.basename(rep).replace('.fa', '')}_{i}\n{species_object.seq}\n")
#     species_file.close()
#     prev_rep = rep
# reference_file.close()


# annotate
def run_prokka(input_fasta):
    shell_command(prokka_cmd.format(os.path.join(os.path.dirname(input_fasta), prokka_version), input_fasta))


def run_eggnog(input_fasta):
    eggnog_dir = os.path.join(os.path.dirname(os.path.dirname(input_fasta)), eggnog_version)
    os.mkdir(eggnog_dir)
    shell_command(eggnog_cmd.format(input_fasta, os.path.join(eggnog_dir, 'eggnog')))


os.chdir(jobs_dir)
sethandlers(file_dir=jobs_dir)

# run prokka
with qp(jobname='prokka', _mem_def='4G', _trds_def=2, _tryrerun=True, _specific_nodes='plink') as q:
    q.startpermanentrun()
    tkttores = {}

    rep_paths = glob.glob(os.path.join(output_dir, '*', '*.fa'))
    for rep_path in rep_paths:
        tkttores[rep_path] = q.method(run_prokka, [rep_path])

    for k, v in tkttores.items():
        q.waitforresult(v)

# run eggnog
with qp(jobname='eggnog', _mem_def='40G', _trds_def=16, _tryrerun=True, _specific_nodes='plink') as q:
    q.startpermanentrun()
    tkttores = {}

    prokka_paths = glob.glob(os.path.join(output_dir, '*', prokka_version, 'prokka.faa'))
    for prokka_path in prokka_paths:
        egg_dir = os.path.join(os.path.dirname(os.path.dirname(prokka_path)), eggnog_version)
        if os.path.exists(os.path.join(egg_dir, 'eggnog.emapper.annotations')):
            pass
        else:
            if os.path.exists(egg_dir):
                os.removedirs(egg_dir)
            tkttores[prokka_path] = q.method(run_eggnog, [prokka_path])

    for k, v in tkttores.items():
        q.waitforresult(v)

# create csvs
csv_path = os.path.join(output_dir, f'{os.path.basename(output_dir)}_{str(datetime.date.today()).replace("-", "_")}_{"{}"}.csv')
prokka_path = os.path.join(output_dir, '*', prokka_version, 'prokka.gff')
eggnog_path = os.path.join(output_dir, '*', eggnog_version, 'eggnog.emapper.annotations')

# create a prokka csv
prokka.OUTPUT_COLUMNS = ['GeneID', 'SGB', 'contig', 'start_pos', 'end_pos', 'strand', 'feature', 'gene', 'product', 'ID']
prokka.gff_to_csv(input_file=glob.glob(prokka_path), output_file=csv_path.format('prokka'))
prokka_df = pd.read_csv(csv_path.format('prokka'))

# create an eggnog csv
shell_command(f"head -5 {glob.glob(eggnog_path)[0]} | tail -1 > {csv_path.format('eggnog')}")  # column names
shell_command(f"grep -v '#' {eggnog_path} >> {csv_path.format('eggnog')}")  # files content
eggnog_df = pd.read_csv(csv_path.format('eggnog'), sep='\t')
eggnog_df['#query_name'] = eggnog_df['#query_name'].str.split(':').str[1]
eggnog_df = eggnog_df.set_index('#query_name')
eggnog_df.to_csv(csv_path.format('eggnog'))

# combine the eggnog and prokka csvs
PE = prokka_df.join(eggnog_df, on='ID')
PE.to_csv(csv_path.format('prokka_eggnog'), index=False)

# zip the run folders
os.chdir(output_dir)
for folder in glob.glob(os.path.join(output_dir, '*')):
    if os.path.isdir(folder) and 'jobs' not in folder:
        folder = os.path.basename(folder)
        shell_command(f"zip {folder}.zip -rm {folder}")
