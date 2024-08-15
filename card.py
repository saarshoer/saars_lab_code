import os
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.GutMBLoader import GutMBLoader

study_ids = ['NICU']
output_dir = '/net/mraid20/export/jasmine/card/NICU'
jobs_dir = os.path.join(output_dir, 'jobs')
lookout_paths = ['/net/mraid20/export/genie/LabData/Data/MBPipeline/NebNext_new/tmp2/UNT/*sample*fastq*',
                 '/net/mraid20/export/genie/LabData/Data/MBPipeline/NovaSeq/tmp2/UNT/*sample*fastq*']

# /net/mraid20/export/genie/LabData/Data/MBPipeline/PipelineRerun/tmp2/HGF
# /net/mraid20/export/genie/LabData/Data/MBPipeline/NebNext_updated/tmp2/UNT

threads = 16
rgi_exe = '/home/saarsh/Develop/Git/rgi/rgi'


def do(f, s):
    os.chdir(output_dir)
    # TODO: make sure we do not want other arguments as well
    os.system(f'python {rgi_exe} bwt -1 {f} -n {threads} -o {s} --clean')
    # TODO: add deletion of unnecessary files


if __name__ == '__main__':

    # queue
    os.makedirs(jobs_dir, exist_ok=True)
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    samples = GutMBLoader().get_data('segal_species', study_ids=study_ids).df_metadata.index.tolist()

    with qp(jobname='rgi', _mem_def='20G', _trds_def=threads,
                _delete_csh_withnoerr=False, _tryrerun=False, _suppress_handlers_warning=True) as q:
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
                    tkttores[sample] = q.method(do, (files[0], sample))
                i = i + 1
                # if i == 2:
                #     break
        print('finished sending jobs')

        print(f'{i}/{len(samples)} ({i / len(samples) * 100}) found')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
