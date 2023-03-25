import os
import numpy as np
import pandas as pd

from PIL import Image
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

from analysis import Study
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE, MAF_1_VALUE

rcParams['font.size'] = 12
rcParams['figure.figsize'] = (
    rcParams['figure.figsize'][0] * 1.6, rcParams['figure.figsize'][1] * 1.6)  # figure size in inches
rcParams['savefig.dpi'] = 200  # figure dots per inch or 'figure'
rcParams['savefig.format'] = 'png'  # {png, ps, pdf, svg}
rcParams['savefig.bbox'] = 'tight'

s = '10K'#'Lifeline_Stool'#'Lifeline_deep'#
run_dir = os.path.join('/net/mraid08/export/jafar/Microbiome/Analyses/saar/strains', s)
strain_var = '/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/mb_sns_strain_variability_R3.h5'
maf_file = os.path.join( f'/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/MAFS/{s}', 'mb_snp_maf_{}_R1_S500.h5')


def do(species):

    os.chdir(run_dir)

    abund = pd.read_pickle(os.path.join('data_frames', 'abundance.df'))

    study = Study(study=s, controls=None, colors=None,
                  alpha=0.05, detection_threshold=0.0001, dissimilarity_threshold=1/20000,
                  base_directory=run_dir)
    study.objs['test'] = study.Object(obj_type='test', df=None, columns=None)

    study.objs['test'].df = pd.read_hdf(os.path.join('data_frames', 'mb_dists.h5'), key=species)
    study.objs['test'].df = study.objs['test'].df.xs('baseline', level='time_point1').\
                                                  xs('baseline', level='time_point2')
    study.objs['test'].df = study.objs['test'].df.join(abund, on=['SampleName1', 'Species'], how='inner').join(
                                                       abund, on=['SampleName2', 'Species'], how='inner',
                                                       lsuffix='1', rsuffix='2')
    study.objs['test'].df = study.objs['test'].df.set_index(['abundance1', 'abundance2'], append=True)

    # dissimilarity
    samples = study.fig_snp_heatmap(obj=study.objs['test'], maximal_filling=0.1, minimal_samples=100,###
                                    method='average',
                                    species=None, cmap='autumn', log_colors=False,
                                    annotations=['abundance'], annotations_cmaps={'abundance': 'binary'},
                                    add_hist=True, add_pairs=False,
                                    strain_var=strain_var)

 #    samples = ['22001804503980_v2_fullrun',
 # '22001804503091_v2_fullrun',
 # '22001804500461_v2_fullrun',
 # '22001804500735_v2_fullrun',
 # '22002003580629_v2_fullrun',
 # '22210480702062_v2_fullrun',
 # '22001804580829_v2_fullrun',
 # '22210580701431_v2_fullrun',
 # '22220180200785_v2_fullrun',
 # '22001804503445_v2_fullrun',
 # '22002003583015_v2_fullrun',
 # '22001804502865_v2_fullrun',
 # '22002003562349_v2_fullrun',
 # '22220180203532_v2_fullrun',
 # '22001804500099_v2_fullrun',
 # '22002003561505_v2_fullrun',
 # '22002003581245_v2_fullrun',
 # '22210580701094_v2_fullrun',
 # '22001804503053_v2_fullrun',
 # '22001804503651_v2_fullrun',
 # '22001804500815_v2_fullrun',
 # '22002010531862_v2_fullrun',
 # '22002010533935_v2_fullrun',
 # '22001804504724_v2_fullrun',
 # '22002003580699_v2_fullrun',
 # '22002003561823_v2_fullrun',
 # '22001804501642_v2_fullrun',
 # '22002003562903_v2_fullrun',
 # '22002003580623_v2_fullrun',
 # '22002010560968_v2_fullrun',
 # '22001804500444_v2_fullrun',
 # '22210580702035_v2_fullrun',
 # '22002003580059_v2_fullrun',
 # '22002003580317_v2_fullrun',
 # '22001804500669_v2_fullrun',
 # '22002003580440_v2_fullrun',
 # '22001804500800_v2_fullrun',
 # '22001804503811_v2_fullrun',
 # '22210580704162_v2_fullrun',
 # '22002003562799_v2_fullrun',
 # '22220180204634_v2_fullrun',
 # '22001804500365_v2_fullrun',
 # '22001804504647_v2_fullrun',
 # '22001804500075_v2_fullrun',
 # '22001804503696_v2_fullrun',
 # '22001804504273_v2_fullrun',
 # '22001804504752_v2_fullrun',
 # '22001804501108_v2_fullrun',
 # '22210580702495_v2_fullrun',
 # '22002003581835_v2_fullrun',
 # '22002010562709_v2_fullrun',
 # '22210580701834_v2_fullrun',
 # '22001804500739_v2_fullrun',
 # '22001804503572_v2_fullrun',
 # '22001804500823_v2_fullrun',
 # '22002003581556_v2_fullrun',
 # '22001804502831_v2_fullrun',
 # '22002010560187_v2_fullrun',
 # '22210480701354_v2_fullrun',
 # '22210480700116_v2_fullrun',
 # '22001804500305_v2_fullrun',
 # '22210680400632_v2_fullrun',
 # '22220180203326_v2_fullrun',
 # '22002010532141_v2_fullrun',
 # '22001804504588_v2_fullrun',
 # '22001804500795_v2_fullrun',
 # '22002003561258_v2_fullrun',
 # '22001804503863_v2_fullrun',
 # '22001804583862_v2_fullrun',
 # '22002003580256_v2_fullrun',
 # '22001804504399_v2_fullrun',
 # '22002003581133_v2_fullrun',
 # '22001804503173_v2_fullrun',
 # '22002003580545_v2_fullrun',
 # '22001804504654_v2_fullrun',
 # '22001804500167_v2_fullrun',
 # '22001804503137_v2_fullrun',
 # '22001804580850_v2_fullrun',
 # '22001804500559_v2_fullrun',
 # '22002003581827_v2_fullrun',
 # '22220180202046_v2_fullrun',
 # '22002003580614_v2_fullrun',
 # '22001804500912_v2_fullrun',
 # '22001804503152_v2_fullrun',
 # '22210480702360_v2_fullrun',
 # '22001804504583_v2_fullrun',
 # '22001804504523_v2_fullrun',
 # '22210580702723_v2_fullrun',
 # '22002010532219_v2_fullrun',
 # '22210580703234_v2_fullrun',
 # '22002003581027_v2_fullrun',
 # '22001804583891_v2_fullrun',
 # '22001804504075_v2_fullrun',
 # '22210480704622_v2_fullrun',
 # '22002003582164_v2_fullrun',
 # '22002003561203_v2_fullrun',
 # '22002010562110_v2_fullrun',
 # '22001804504352_v2_fullrun',
 # '22220180203926_v2_fullrun',
 # '22001804504469_v2_fullrun',
 # '22002003562642_v2_fullrun',
 # '22002010532742_v2_fullrun',
 # '22220180203363_v2_fullrun',
 # '22002003581078_v2_fullrun',
 # '22001804503418_v2_fullrun',
 # '22002010533404_v2_fullrun',
 # '22220180203466_v2_fullrun',
 # '22220180203809_v2_fullrun',
 # '22002003583017_v2_fullrun',
 # '22001804502946_v2_fullrun',
 # '22002003580786_v2_fullrun',
 # '22001804504743_v2_fullrun',
 # '22002010532361_v2_fullrun',
 # '22002010530414_v2_fullrun',
 # '22002010531350_v2_fullrun',
 # '22002010532584_v2_fullrun',
 # '22001804503168_v2_fullrun',
 # '22001804503716_v2_fullrun',
 # '22002003580429_v2_fullrun',
 # '22210580702648_v2_fullrun',
 # '22001804583856_v2_fullrun',
 # '22001804504524_v2_fullrun',
 # '22001804504637_v2_fullrun',
 # '22001804580882_v2_fullrun',
 # '22210580702994_v2_fullrun',
 # '22001804500603_v2_fullrun',
 # '22210480703393_v2_fullrun',
 # '22001804504686_v2_fullrun',
 # '22210580700528_v2_fullrun',
 # '22001804504786_v2_fullrun',
 # '22001804502138_v2_fullrun',
 # '22002010531470_v2_fullrun']

    if samples is not None:

        # shared_pos
        drop = list(set(study.objs['test'].df.index.names) - set(['SampleName1', 'SampleName2']))
        plt.figure(figsize=(8.9, 8.9))
        ax = sns.heatmap(study.objs['test'].df['shared_pos'].reset_index(drop, drop=True).
                        unstack('SampleName2').loc[samples, samples],
                        cmap='autumn', cbar_kws={'label': 'shared positions', 'orientation': 'horizontal',
                                                 'ticks': None, 'shrink': 0.55, 'pad': 0.055, 'aspect': 40},
                        xticklabels=False, yticklabels=False)
        plt.xlabel('SampleName1')
        plt.ylabel('SampleName2')
        ax.yaxis.set_label_position('right')
        plt.savefig(os.path.join('figs', 'test', f'{species}_shared_positions'), bbox_inches='tight')
        plt.close()

        # MAF
        position_threshold = 0.95
        df = None

        with pd.HDFStore(maf_file.format(species), 'r') as hdf:
            for key in hdf.keys():
                contig_part_df = hdf[key].fillna(MAF_MISSING_VALUE)
                contig_part_df = contig_part_df[contig_part_df.index.get_level_values('SampleName').isin(samples)]
                if type(contig_part_df) == pd.DataFrame:
                    contig_part_df = contig_part_df.stack()
                contig_part_df.index = contig_part_df.index.reorder_levels(['Position', 'SampleName'])
                variable_positions = ((contig_part_df.unstack('Position') == MAF_1_VALUE).sum() <=
                                      len(samples)*position_threshold).replace(False, np.nan).dropna().index.tolist()
                missing_positions = ((contig_part_df.unstack('Position') == MAF_MISSING_VALUE).sum() <=
                                     len(samples)*position_threshold).replace(False, np.nan).dropna().index.tolist()
                contig_part_df = contig_part_df[contig_part_df.index.get_level_values('Position').isin(
                    variable_positions+missing_positions)]
                contig_part_df = contig_part_df.to_frame().assign(contig=int(key.split('_')[1]))

                df = pd.concat([df, contig_part_df]) if df is not None else contig_part_df.copy()

            df = df.set_index('contig', append=True).unstack(['contig', 'Position']).replace(MAF_MISSING_VALUE, np.nan)
            df.columns = df.columns.droplevel(0)
            df = df.loc[samples, sorted(df.columns)]

            plt.figure(figsize=(8.9, 8.9))
            ax = sns.heatmap(df, cmap='autumn', cbar_kws={'label': 'MAF', 'orientation': 'horizontal',
                                                         'ticks': None, 'shrink': 0.55, 'pad': 0.055, 'aspect': 40},
                            xticklabels=False, yticklabels=False)
            plt.xlabel('Contig-Position')
            plt.ylabel('SampleName2')
            ax.yaxis.set_label_position('right')
            plt.savefig(os.path.join('figs', 'test', f'{species}_MAF'), bbox_inches='tight')
            plt.close()

        # summary
        images = [Image.open(x) for x in [f'figs/test/{species}.png',  # order matters
                                          f'figs/test/{species}_MAF.png',
                                          f'figs/test/{species}_shared_positions.png']]
        widths, heights = zip(*(i.size for i in images))
        new_im = Image.new(mode='RGB', size=(sum(widths)-400, max(heights)), color='white')
        x_offset = 0
        for i, im in enumerate(images):
            new_im.paste(im, (x_offset - (0 if i == 0 else 400), (0 if i == 0 else 670)))
            x_offset += im.size[0]
        new_im.save(os.path.join('figs', 'test', f'{species}_summary.png'))

    os.chdir(jobs_dir)


if __name__ == '__main__':

    jobs_dir = os.path.join(run_dir, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    with pd.HDFStore(os.path.join(run_dir, 'data_frames', 'mb_dists.h5'), 'r') as hdf:
        keys = list(hdf.keys())
    del hdf

    with qp(jobname='test', _tryrerun=True, _mem_def='20G') as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        for key in keys:#['/Rep_1551']:#
            tkttores[key] = q.method(do, (key[1:],), _job_name=f't{key.split("_")[-1]}')
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

    print('done')
