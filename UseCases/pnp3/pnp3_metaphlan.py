import os
import pandas as pd
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

input_files = pd.read_csv(
    '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/PNP3/paper/microbiome_samples_location.csv').iloc[:, 0].tolist()
output_dir = '/net/mraid20/ifs/wisdom/segal_lab/jafar/Microbiome/Analyses/saar/PNP3/paper/microbiome_metaphlan/'


def do(in_file, out_file):
    mpa_command = '/usr/wisdom/python/bin/metaphlan --input_type fastq --bowtie2db /net/mraid20/export/genie/Data/Databases/mpa_2023 --no_map -t rel_ab'
    mpa_command = ' '.join([mpa_command, in_file, out_file])
    os.system(mpa_command)


if __name__ == '__main__':
    jobs_dir = os.path.join(os.path.dirname(output_dir), 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    with qp(jobname='metaphlan', _tryrerun=False, _delete_csh_withnoerr=True, _mem_def='20G', _trds_def=4, max_r=100) as q:
        # _tryrerun=False otherwise some metaphlan info is interpreted as error

        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for in_f in input_files:
            out_f = os.path.join(output_dir, os.path.basename(in_f).replace('.fastq.gz', '.tsv'))
            if not os.path.exists(out_f):
                tkttores[out_f] = q.method(do, [in_f, out_f])
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')
