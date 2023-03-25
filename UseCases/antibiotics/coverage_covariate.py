import os
import glob
import pandas as pd
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

study = '10K'
run_type = 'within'  # always

run_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/{run_type}'
counts_dir = f'/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/CountsS/'


def do(f):
    species = f'Rep_{f.split("_")[-1].split(".")[0]}'
    samples = pd.read_pickle(os.path.join(os.path.dirname(run_dir), 'data_frames', 'snps.df')).index.tolist()
    positions = pd.read_hdf(f).reset_index().set_index('Contig')['Position']
    positions.index = ['/' + '_'.join(contig_part.split('_')[:2]) for contig_part in positions.index]

    positions_len = 0
    positions_sum = pd.Series(index=samples).fillna(0)

    with pd.HDFStore(os.path.join(os.path.dirname(run_dir), 'coverage_covariate', f'mb_snp_coverage_{species}.h5'),
                     complevel=9) as coverage:

        with pd.HDFStore(os.path.join(counts_dir, f'mb_snp_counts_{species}.h5'), 'r') as counts:

            for contig_part in counts.keys():
                contig = '_'.join(contig_part.split('_')[:2])
                if contig not in positions.index:
                    continue
                contig_positions = positions.loc[contig].tolist()
                if type(contig_positions) is int:
                    contig_positions = [contig_positions]

                contig_counts = counts[contig_part].sum(axis=1)
                contig_counts = contig_counts.loc[contig_counts.index.get_level_values('SampleName').isin(samples)]
                contig_counts = contig_counts.loc[contig_counts.index.get_level_values('Position').isin(contig_positions)]
                contig_counts = contig_counts.unstack('Position')

                if contig_counts.empty:
                    continue

                try:
                    coverage[contig] = coverage[contig].join(contig_counts, how='outer')
                except ValueError:
                    print(contig)
                    print(f'problematic positions {contig_counts.columns[contig_counts.columns.isin(coverage[contig].columns)]}')
                    contig_counts = contig_counts.loc[:, ~contig_counts.columns.isin(coverage[contig].columns)]
                    coverage[contig] = coverage[contig].join(contig_counts, how='outer')
                except KeyError:
                    coverage[contig] = contig_counts

                positions_len = positions_len + contig_counts.shape[1]
                positions_sum = positions_sum.add(contig_counts.sum(axis=1).fillna(0), fill_value=0)

        positions_mean = positions_sum/positions_len
        with pd.HDFStore(os.path.join(os.path.dirname(run_dir), 'coverage_covariate', f'mb_snp_covariate_{species}.h5'),
                         complevel=9) as covariate:

            for contig in coverage.keys():
                positions_mean_contig_samples = positions_mean.loc[coverage[contig].index]
                covariate[contig] = coverage[contig].apply(lambda col: col / positions_mean_contig_samples)

    if positions_len != positions.shape[0]:
        print(f'got {positions_len:,}, expected {positions.shape[0]:,}')


if __name__ == '__main__':

    jobs_dir = os.path.join(run_dir, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(os.path.join(os.path.dirname(run_dir), 'coverage_covariate'))

    with qp(jobname='covcov', _tryrerun=True, _delete_csh_withnoerr=False, _mem_def='30G') as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        tested = glob.glob(os.path.join(run_dir, 'naive_hdfs', 'mb_gwas_Rep_*_Rep_*.h5'))
        for file in tested:
            tkttores[file] = q.method(do, [file], _job_name=f'c{file.split("_")[-1].split(".")[0]}')
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
