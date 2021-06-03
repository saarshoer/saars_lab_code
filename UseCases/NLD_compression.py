import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabMBPipeline.ProcessUtils import _shell_command

folders = [
    '/net/mraid08/export/genie/Analyses/Dists/Data/NLD1/older/bams/',
    '/net/mraid08/export/genie/Analyses/Dists/Data/NLD1/older/fastqs/',
    '/net/mraid08/export/genie/Analyses/Dists/Data/NLD2/bams/',
    '/net/mraid08/export/genie/Analyses/Dists/Data/NLD2/fastqs/',
    '/net/mraid08/export/genie/Analyses/Dists/Data/NLD2/fastqs/split/',
]

# '/net/mraid08/export/genie/Analyses/Dists/Data/NLD1/older/gzip/',
# '/net/mraid08/export/genie/Analyses/Dists/Data/NLD2/gzip2/'

jobs_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/NLDcopmJobs'


def func(folder):
    cmd = f'rm -Rf {folder} &'
    print(cmd)
    _shell_command(cmd)

    # files = glob.glob(os.path.join(folder, '*'))
    # len_files = len(files)
    # for i_file, file in enumerate(files):
    #     if not os.path.isdir(file):
    #         print(f'file {i_file + 1}/{len_files}')
    #             _shell_command('gzip -9 ' + file)


# queue
os.chdir(jobs_dir)
sethandlers(file_dir=jobs_dir)

with qp(jobname='NLDcomp', _mem_def='10G', _tryrerun=False) as q:
    q.startpermanentrun()
    tkttores = {}

    for i_folder, folder in enumerate(folders):
        tkttores[i_folder] = q.method(func, [folder])

    for k, v in tkttores.items():
        q.waitforresult(v)
