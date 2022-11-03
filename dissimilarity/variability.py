import os
import glob
import numpy as np
import pandas as pd
from LabQueue.qp import qp
from LabUtils.addloglevels import sethandlers


def func(species_fname, output_fname):

    variable_pos = total_pos = pd.Series(dtype=np.float64)

    with pd.HDFStore(species_fname, 'r') as hdf:

        for contig_part in hdf.keys():
            contig_df = hdf[contig_part]
            contig_df = contig_df[contig_df.sum(axis=1) >= min_reads]

            variable_pos = variable_pos.add(
                pd.DataFrame(contig_df.groupby('SampleName').apply(lambda s: ((s != 0).sum(axis=1) > 1).sum())),
                fill_value=0)
            total_pos = total_pos.add(
                pd.DataFrame(contig_df.groupby('SampleName').apply(lambda s: s.shape[0])),
                fill_value=0)

    df = variable_pos.div(total_pos).rename(
         columns={0: f'{species3letters}_{species_fname.split("_")[-1].split(".")[0]}'})
    df.to_hdf(output_fname, key='snps')


if __name__ == '__main__':

    # PNP3
    # counts_dir = '/net/mraid08/export/genie/LabData/Data/MBPipeline/PNP3_rerun_segata/MBSNP/CountsB/'
    # species3letters = 'SGB'

    # Everything else
    counts_dir = '/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/CountsS/'
    species3letters = 'Rep'

    counts_fname = f'mb_snp_counts_{species3letters}_*.h5'
    summary_fname = 'mb_snb_strain_variability_R{}.h5'

    min_reads = 3

    jobs_dir = os.path.join(os.path.dirname(os.path.dirname(counts_dir)), 'jobs_strains')

    sethandlers(file_dir=jobs_dir)
    os.chdir(jobs_dir)

    with qp(jobname='strains', _delete_csh_withnoerr=True, q=['himem7.q'], _mem_def='10G',
            delay_batch=10, delay_sec=-1) as q:

        q.startpermanentrun()
        tkttores = {}

        species_files = glob.glob(os.path.join(counts_dir, counts_fname))
        for i, species_fname in enumerate(species_files):
            output_fname = os.path.join(jobs_dir, os.path.basename(species_fname))
            # if not os.path.exists(output_fname):
            #     print(output_fname)
            tkttores[output_fname] = q.method(func, [species_fname, output_fname])
            # if i == 3:
            #     break

        summary_df = pd.DataFrame()
        for output_fname, v in tkttores.items():
            q.waitforresult(v)

        # for species_fname in species_files:
        #     output_fname = os.path.join(jobs_dir, os.path.basename(species_fname))

            summary_df = summary_df.join(pd.read_hdf(output_fname), how='outer')

        summary_df.to_hdf(os.path.join(os.path.dirname(os.path.dirname(counts_dir)), summary_fname.format(min_reads)),
                          key='snps')

    #     for output_fname, v in tkttores.items():
    #         os.remove(output_fname)
    #
    print('done')
