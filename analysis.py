# Imports
# files
import os
import pickle
import tarfile
from glob import glob
from ftplib import FTP

# data
import copy
import random
import numpy as np
import pandas as pd
from itertools import combinations

# statistics
from skbio import DistanceMatrix
from scipy.spatial import distance_matrix
from mne.stats.multi_comp import bonferroni_correction, fdr_correction
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_rel, ttest_1samp, binom_test, \
    pearsonr, spearmanr, kendalltau
from statannot import add_stat_annotation

# models
from statsmodels.api import Logit
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier, XGBRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, accuracy_score
from sklearn.model_selection import RepeatedKFold, train_test_split, GridSearchCV

# dimensionality reduction
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# clustering
from scipy.signal import find_peaks
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster

# plots
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams, cm, colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# lab
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders import GutMBLoader, OralMBLoader, BloodTestsLoader, BodyMeasuresLoader, \
    DietLoggingLoader

segata_df = pd.read_csv('/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Segata/'
                        'SupplementaryTable8-SGBsDescription.csv')

segal_df = pd.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/saar/GutMBLoader.xlsx').set_index('sBin')


# Main class
class Study:

    # Initialization
    def __init__(self,

                 # Default parameters
                 study=None, controls=None, colors=None,
                 alpha=0.05, detection_threshold=0.0001, dissimilarity_threshold=1/20000,
                 base_directory=os.getcwd(), fig_fontsize=12, fig_sizefactor=1.5, fig_dpi=200, fig_format='png'):

        # Parameters
        self.params = _Parameters(study=study, controls=controls, colors=colors,
                                  alpha=alpha,
                                  detection_threshold=detection_threshold,
                                  dissimilarity_threshold=dissimilarity_threshold)

        # Directories
        self.dirs = _Directories(base_directory)

        # Objects
        self.objs = {}

        # Figures
        # full list of parameters in https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
        rcParams['font.size'] = fig_fontsize
        rcParams['figure.figsize'] = (rcParams['figure.figsize'][0] * fig_sizefactor,
                                      rcParams['figure.figsize'][1] * fig_sizefactor)  # figure size in inches
        rcParams['savefig.dpi'] = fig_dpi  # figure dots per inch or 'figure'
        rcParams['savefig.format'] = fig_format  # {png, ps, pdf, svg}

        # Jobs
        sethandlers(file_dir=self.dirs.jobs)

    class Object:

        # Initialization
        def __init__(self, obj_type=None, df=None, columns=None):
            self.type = obj_type
            self.df = df
            self.columns = columns

    # Functions

    # data generation
    def load_df(self, data,
                columns_from_metadata=None,
                columns_from_file=None, file=None, file_index=None):
        """
        Load a data frame and add columns from external sources

        :param data: (str or pd.DataFrame) string of data type ('gutMB', 'oralMB', 'blood', 'body', 'cgm', 'diet')
         or a data frame
        :param columns_from_metadata: (str, list or dict) to add from the LabDataLoader meta data data frame
        :param columns_from_file: (str, list or dict) to add from a file
        :param file: (str or pd.DataFrame) from which to add columns ('xlsx' or 'csv')
        :param file_index: (str, list or dict) index in file to join by

        :return: (pd.DataFrame) df
        """

        def fig_abundance_reads(input_metadata_df):
            columns = ['Raw', 'PostQC', 'PostHGF', 'Unaligned']  # column names without 'RC'

            # remove 'RC' from column name because y-axis label is read length
            input_metadata_df.columns = [col_.replace('RC', '') for col_ in input_metadata_df.columns]

            plt.figure()
            sns.boxplot(data=input_metadata_df[columns])
            title = '{} reads count'.format(data.replace('MB', ''))
            plt.title(title)
            plt.ylabel('reads count')

            plt.savefig(os.path.join(self.dirs.figs, title))

            # resizing over PNP3 gut outliers
            if self.params.study == 'PNP3' and data == 'gutMB':
                plt.ylim([0, 4*(10**7)])
                plt.savefig(os.path.join(self.dirs.figs, '{} - limited to 40M'.format(title)))

        # default parameters to always act on
        default_indices_order = ['group', 'person', 'time_point', 'time', 'sample']
        default_name_changes = {'SampleName': 'sample', 'RegistrationCode': 'person', 'Date': 'time', 'alloc': 'group'}

        # loading from Lab Data
        if self.params.study == 'PNP3':
            self.params.study = ['PNP3', 'T2D']
            reg_ids = '/net/mraid08/export/genie/LabData/Data/StudySpecificData/PNP3/Combined_RandomizationResults_final.xlsx'
            reg_ids = pd.read_excel(reg_ids, engine='openpyxl')
            reg_ids = reg_ids.loc[reg_ids['IsActive'] == 1, 'RegistrationCode'].astype(int).astype(str).tolist()
        else:
            reg_ids = None

        if type(data) == str:  # TODO: can be modified to get a loader and its params as arguments
            if data == 'gutMB':
                lab_data = GutMBLoader.GutMBLoader().get_data('segata_species', study_ids=self.params.study, reg_ids=reg_ids)
                fig_abundance_reads(lab_data.df_metadata)
            elif data == 'oralMB':
                lab_data = OralMBLoader.OralMBLoader().get_data('segata_species', study_ids=self.params.study, reg_ids=reg_ids)
                fig_abundance_reads(lab_data.df_metadata)
            elif data == 'blood':
                lab_data = BloodTestsLoader.BloodTestsLoader().get_data(study_ids=self.params.study, reg_ids=reg_ids)
            elif data == 'body':
                lab_data = BodyMeasuresLoader.BodyMeasuresLoader().get_data(study_ids=self.params.study, reg_ids=reg_ids)
            elif data == 'cgm':
                pass
                # lab_data = CGMLoader.CGMLoader().get_data(study_ids=self.params.study)
            elif data == 'diet':
                loader = DietLoggingLoader.DietLoggingLoader()
                lab_data = loader.get_data(study_ids=self.params.study, reg_ids=reg_ids)
                lab_data.df = loader.add_nutrients(lab_data.df,
                                                   nutrient_list=['protein_g', 'totallipid_g', 'carbohydrate_g'])
            else:
                raise Exception('data string is not valid')

            if reg_ids:
                self.params.study = 'PNP3'

            # retrieving
            df = lab_data.df
            metadata_df = lab_data.df_metadata

        # loading from data frame
        elif type(data) == pd.DataFrame:
            df = data
            metadata_df = None  # just to prevent a warning later on

        # unrecognizable data type
        else:
            raise Exception('data type is not valid')

        added_columns = []

        # add to df columns from metadata
        if columns_from_metadata is not None and metadata_df is not None:

            # allow the user to enter columns_from_metadata as a list and as a
            # dictionary with column name replacement or str
            if type(columns_from_metadata) is dict:
                metadata_df = metadata_df.rename(columns=columns_from_metadata)
                columns_from_metadata = list(columns_from_metadata.values())
            elif type(columns_from_metadata) is str:
                columns_from_metadata = [columns_from_metadata]
            elif type(columns_from_metadata) is not list:
                raise Exception('columns_from_metadata type is not valid')

            added_columns = added_columns + columns_from_metadata

            # add meta data columns to each entry, if is missing leave out
            for col in columns_from_metadata:
                df = df.join(metadata_df[col].dropna(), on=metadata_df.index.names, how='inner')

        # add to df columns from file
        if columns_from_file is not None and file is not None:

            if type(file) is str:
                # read file to data frame according to the file's type
                file_extension = os.path.splitext(file)[1]
                if file_extension == '.xlsx':
                    file_df = pd.read_excel(file, engine='openpyxl')
                elif file_extension == '.csv':
                    file_df = pd.read_csv(file)
                else:
                    raise Exception('file extension is not valid')

            # loading from data frame
            elif type(file) is pd.DataFrame:
                file_df = file.reset_index()

            # unrecognizable file type
            else:
                raise Exception('file type is not valid')

            # fix known problems
            if 'RegistrationCode' in file_df.columns:
                file_df['RegistrationCode'] = file_df['RegistrationCode'].astype(str).str.split('.').str[0]

            if 'Timestamp' in file_df.columns:
                file_df['Timestamp'] = pd.to_datetime(file_df['Timestamp'], dayfirst=True, errors='coerce')

            # allow the user to enter columns_from_file as a list and as a
            # dictionary with column name replacement or str
            if type(columns_from_file) is dict:
                file_df = file_df.rename(columns=columns_from_file)
                columns_from_file = list(columns_from_file.values())
            elif type(columns_from_file) is str:
                columns_from_file = [columns_from_file]
            elif type(columns_from_file) is not list:
                raise Exception('columns_from_file type is not valid')

            added_columns = added_columns + columns_from_file

            # set the file_df index to be whatever index you want to match on with df
            if type(file_index) is dict:
                file_df = file_df.rename(columns=file_index)
                file_index = list(file_index.values())
            elif type(file_index) is str:
                file_index = [file_index]
            elif type(file_index) is not list:
                raise Exception('file_index type is not valid')

            # remove duplicate rows only by the relevant index/columns for cases where tests are in rows
            file_df = file_df.drop_duplicates(subset=file_index + columns_from_file)

            # set index
            file_df = file_df.set_index(file_index)

            # add file columns to each entry, if is missing leave out
            for col in columns_from_file:
                df = df.join(file_df[col].dropna(), on=file_df.index.names, how='inner')

            # TODO: metadata and file handling can be moved to functions

        # index manipulations
        df = df.set_index(added_columns, append=True)
        df = df.rename_axis(index=default_name_changes)

        indices_order = list(df.index.names)
        indices_order.sort(key=lambda i: (default_indices_order + list(df.index.names)).index(i))
        df = df.reorder_levels(indices_order)
        df = df.sort_values(by=df.index.names)

        return df

    def prep_df(self, obj, time_points=None, mean_time_point=False, n_entries=None, indices_dict=None,
                log=10, clip_std_quantile=False, robust_standardization=False, clip_std=False, fill_na=False):
        """
        Handel the time_point, replace indices values, handel missing values (if abundance) and optionally log10

        :param obj:
        :param time_points: (list) to include
        :param mean_time_point: (bool) if to mean between duplicate entries of the same time_point
        :param n_entries: (int) amount of entries to accept per person
        :param indices_dict: (dict of dicts) outer key for index name
        inner key for index value and value after replacement
        :param log: (int or str or bool) [2, 10, 'e' or False] whether to log
        :param clip_std_quantile: (float) how many quantile to trim from each side
        :param robust_standardization: (bool) whether to robust standardize the data (using median)
        :param clip_std: (int) how many standard deviations to clip from each side
        :param fill_na: (bool) whether to fill missing values with column's minimal value (True) or other given value


        :return: (pd.DataFrame) df
        """

        def fig_abundance_distribution(abundance_df):

            df1 = pd.DataFrame(abundance_df.stack()).reset_index()
            df2 = pd.DataFrame(np.log10(abundance_df).stack()).reset_index()
            df1['log'] = False
            df2['log'] = True

            # if there is no actual control time_point we have a problem
            if self.params.controls['time_point'] != 'fake_time_point':
                df3 = pd.DataFrame(get_delta_df(abundance_df, self.params.controls['time_point'])
                                   .stack()).reset_index()
                df4 = pd.DataFrame(get_delta_df(np.log10(abundance_df), self.params.controls['time_point'])
                                   .stack()).reset_index()
                df3['log'] = False
                df4['log'] = True
                full_df = pd.concat([df1, df2, df3, df4], sort=False).rename(columns={0: 'abundance'})
            else:
                full_df = pd.concat([df1, df2], sort=False).rename(columns={0: 'abundance'})
                full_df['time_point'] = 'fake_time_point'

            g = sns.FacetGrid(full_df, row='log', col='time_point', hue='group', palette=self.params.colors,
                              sharex=False, sharey=False, margin_titles=True)
            g = g.map(sns.distplot, 'abundance', hist=True, kde=False)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            title = '{} distribution'.format(obj.type)
            plt.suptitle(title)
            plt.subplots_adjust(top=0.9)

            plt.savefig(os.path.join(self.dirs.figs, title))

        def fig_species_distribution(abundance_df):
            n_species_per_sample = pd.DataFrame(pd.DataFrame(
                (abundance_df > self.params.detection_threshold).sum(axis=1))
                                                .stack()).reset_index().rename(columns={0: f'n_{obj.columns}'})

            if self.params.controls['time_point'] == 'fake_time_point':
                n_species_per_sample['time_point'] = 'fake_time_point'

            # TODO: replace repeating code for 'group' and for 'time_point' by major and minor
            # hue 'group'
            g = sns.FacetGrid(n_species_per_sample, col='time_point', hue='group', palette=self.params.colors,
                              sharex=True, sharey=True, margin_titles=True)
            g = g.map(sns.distplot, f'n_{obj.columns}', hist=True, kde=False, bins=10)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            g.set_ylabels('n_samples')
            title = '{} distribution colored by group'.format(obj.type.replace('abundance', '')).replace('  ', ' ')
            plt.suptitle(title)
            plt.subplots_adjust(top=0.8)

            plt.savefig(os.path.join(self.dirs.figs, title))

            # hue 'time_point'
            g = sns.FacetGrid(n_species_per_sample, col='group', hue='time_point', palette=self.params.colors,
                              sharex=True, sharey=True, margin_titles=True)
            g = g.map(sns.distplot, f'n_{obj.columns}', hist=True, kde=False, bins=10)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            g.set_ylabels('n_samples')
            title = '{} distribution colored by time_point'.format(obj.type.replace('abundance', '')).replace('  ', ' ')
            plt.suptitle(title)
            plt.subplots_adjust(top=0.8)

            plt.savefig(os.path.join(self.dirs.figs, title))

        def time2time_point(person_abundance_df):

            # add time_point by the time of the points, for example: 03/05, 01/04 will become 1, 0
            person_abundance_df['time_point'] = person_abundance_df.index.get_level_values('time').argsort()

            return person_abundance_df

        def declare_missing_values(person_abundance_df):

            # in places where all column values equal to detection threshold declare them as missing values
            person_abundance_df.loc[:, person_abundance_df.isna().all()] = 'NA'

            return person_abundance_df

        # retrieving
        df = obj.df

        # add time_point by the time of the points, for example: 03/05, 01/04 will become 1, 0
        if 'time_point' not in df.index.names and 'time' in df.index.names:
            df = df.groupby(['person']).apply(time2time_point)
            indices_order = list(df.index.names)
            df = df.set_index('time_point', append=True)
            indices_order.insert(indices_order.index('time') + 1, 'time_point')
            df = df.reorder_levels(indices_order)

        # take only the wanted time_points
        if time_points is not None:
            df = df[df.index.get_level_values('time_point').isin(time_points)]

        # average multiple entries from the same time
        if mean_time_point:
            indices = list(df.index.names)
            indices.remove('time')
            df = df.reset_index().groupby(['person', 'time_point']).mean().reset_index()
            df = df.set_index(indices)
            # obliviating time index because the mean function can not handel datetime values

        # filter out cases where a person has the wrong number of samples
        if n_entries is not None:
            df = df.groupby('person').filter(lambda x: len(x) == n_entries)

        # replace indices values according to dictionary
        if indices_dict is not None:
            for ind in indices_dict.keys():
                if ind in df.index.names:
                    df = df.rename(index=indices_dict[ind], level=ind)

        # filter out empty columns (species/test)
        df = df.dropna(axis=1, how='all')

        if 'abundance' in obj.type:
            # figures
            fig_abundance_distribution(df)
            fig_species_distribution(df)

        # log
        if log:
            log_funcs = {10: np.log10, 2: np.log2, 'e': np.log}
            df = log_funcs[log](df.astype(float))

        # robust standardization
        if robust_standardization:
            # clip just the standard deviation's quantile
            if clip_std_quantile:
                clipped_std = df.clip(lower=df.quantile(clip_std_quantile),
                                      upper=df.quantile(1 - clip_std_quantile), axis=1).std()
            else:
                clipped_std = df.std()
            df = (df - df.median()) / clipped_std

        # clip standard deviation
        if clip_std:
            for col in df.columns:
                col_mean = df[col].mean()
                col_std = df[col].std()
                df[col] = df[[col]].clip(lower=col_mean - clip_std * col_std,
                                         upper=col_mean + clip_std * col_std, axis=1)

        # fill missing values with minimum
        if fill_na:

            if fill_na is True:
                na_value = df.min()
            else:
                na_value = fill_na

            # if all person's values are missing it will not be filled
            if len(np.unique(df.index.get_level_values('time_point'))) > 1:
                df = df.groupby(['person']).apply(declare_missing_values)

            df = df.fillna(na_value).replace('NA', np.nan)

        return df

    # analysis
    def comp_stats(self, obj, test, between, delta=False, minimal_samples=0.1, main_correction_method='FDR',
                   normalize_figure=False, internal_use=False):
        """
        Compute the statistical significance for the difference between elements

        :param obj: (object) object to calculate the statistics on
        :param test: (str) mannwhitneyu, wilcoxon, ttest_ind, ttest_rel, ttest_1samp, binom_test
        :param between: (str) group or time_point
        :param delta: (bool) whether or not to calculate the statistics based on the difference between
        time_point values
        :param minimal_samples: (percentage or int) to have in order to compute
        :param main_correction_method: (str) bonferroni or FDR, main multiple hypothesis correction method
        :param normalize_figure: (bool) whether to normalizes values scale in figure
        :param internal_use: (bool) whether an internal function is running it or not

        :return: None or pd.DataFrame
        """

        # TODO: remove internal use

        # TODO: maybe in case of abundance, remove detection threshold form comparisons
        #  that do not include delta or paired

        def fig_significant_stats(figure_internal, axes_internal, curr_data_df):

            # first run
            if major_e == major_elements[0] and minor_e == minor_elements[0]:

                # create figure object
                if len(major_elements) == 1 and len(minor_elements) == 1:  # in case of a single plot
                    figure_internal, axes_internal = \
                        plt.subplots(figsize=(rcParams['figure.figsize'][0]/2, rcParams['figure.figsize'][1]))
                    axes_internal = [axes_internal]

                else:  # in case of multiple plots
                    figure_internal, axes_internal = \
                        plt.subplots(nrows=min(len(major_elements), len(minor_elements)),
                                     ncols=max(len(major_elements), len(minor_elements)),
                                     sharey=True)

                plt.subplots_adjust(bottom=0.5, left=0.14)
                plt.suptitle('{}\n{} - significant results'.format(obj.type, test))

            ax = axes_internal[major_elements.index(major_e) * len(minor_elements) + minor_elements.index(minor_e)]

            # in case there are some significant results
            curr_data_df = curr_data_df.dropna()
            if curr_data_df.shape[0] != 0:

                # transfer the columns to rows with identifiers
                # reset index behaves differently for single level versus multilevel indices
                # TODO: handle case where input data frame has multiple level columns
                if len(curr_data_df.index.names) == 1 and curr_data_df.index.names != [None]:
                    name2replace = curr_data_df.index.names[0]
                else:
                    name2replace = 'level_0'
                curr_data_df = pd.DataFrame(curr_data_df.stack()).reset_index() \
                    .rename(columns={name2replace: 'index', 0: 'list'})

                # split the lists to a row for each element in them
                curr_data_df = curr_data_df \
                    .merge(curr_data_df['list'].apply(pd.Series), right_index=True, left_index=True) \
                    .melt(id_vars=['index', major, 'list'], value_name='value') \
                    .drop(['list', 'variable'], axis=1) \
                    .dropna()

                # normalize in case measurements are in different scales
                if normalize_figure:
                    curr_data_df['value'] = curr_data_df.groupby('index').transform(
                        lambda x: (x-x.min()) / (x.max()-x.min()))
                    # TODO: fix cases where delta is True

                # clean the bacteria name
                if 'abundance' in obj.type and obj.columns == 'bacteria':
                    curr_data_df['index'] = curr_data_df['index'].apply(
                            lambda row: '{} ({})'.format(row.split('|')[6], row.split('|')[-1]).replace('s__', ''))
                    y_label = 'abundance'
                else:
                    y_label = '{} values'.format(obj.columns)

                # clean the blood test name
                if 'blood' in obj.type:
                    curr_data_df['index'] = curr_data_df['index'].apply(
                        lambda row: row.replace('bt__', ''))

                # clean the nutrient name
                if 'diet' in obj.type:
                    curr_data_df['index'] = curr_data_df['index'].apply(
                        lambda row: row.replace('_g', ''))

                # sort index values so they will be sorted in the plot
                curr_data_df = curr_data_df.sort_values(by='index')

                # plotting
                sns.boxplot(x='index', y='value', hue=major, hue_order=np.unique(curr_data_df[major]),
                            palette=self.params.colors, data=curr_data_df, ax=ax)
                ax.set_ylabel('normalized {}'.format(y_label)) if normalize_figure else ax.set_ylabel(y_label)
                ax.set_xlabel('')
                ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=90)
                ax.legend().set_visible(False)
                if delta:
                    xmin, xmax = ax.get_xlim()
                    ax.hlines(y=0, xmin=xmin, xmax=xmax, colors='lightgray')

            else:
                # texting
                ax.text(x=0.5, y=0.5, s='no significant results', horizontalalignment='center', transform=ax.transAxes)
                ax.axis('off')

            ax.set_title(minor_e, color=self.params.colors[minor_e])

            # last run
            if major_e == major_elements[-1] and minor_e == minor_elements[-1]:

                # legend
                label_list = []
                handle_list = []
                for ax in axes_internal:
                    handles, labels = ax.get_legend_handles_labels()
                    for handle, label in zip(handles, labels):
                        if label not in label_list:
                            label_list.append(label)
                            handle_list.append(handle)
                figure_internal.legend(handle_list, label_list, loc='lower center')

                # saving
                plt.savefig(os.path.join(
                    self.dirs.figs,
                    'stats {}{} {} {}'.format(obj.type, delta_str, test, between)))
                plt.close()
                figure_internal = None
                axes_internal = None

            return figure_internal, axes_internal

        # retrieving the data frames from the objects
        df = obj.df

        # change the values in the data frames to be the change in value between two time points:
        # any time point and the control time point
        if delta:
            df = get_delta_df(df, self.params.controls['time_point'])
            delta_str = ' delta'  # for file names
        else:
            delta_str = ''  # for file names

        # convert the minimal_samples percentage from the argument to a minimal_samples as number of samples
        if minimal_samples < 1:  # meaning minimal_samples is in percentage
            minimal_samples = round(minimal_samples * df.shape[0])

        # the tests can be between groups (control/test1/test2/...) or between time points (0/1/2...)
        between_options = ['group', 'time_point']  # this list has to be length 2
        if between not in between_options:
            raise Exception('between not valid')

        # names of elements
        major = between
        minor = between_options[abs(between_options.index(major) - 1)]  # the index of the not major element

        # list of elements
        major_elements = np.unique(df.index.get_level_values(major)).tolist()
        minor_elements = np.unique(df.index.get_level_values(minor)).tolist()

        # controls values
        major_control = self.params.controls[major]

        # tests without control data
        tests_without_control = ['ttest_1samp', 'binom_test']

        # data frame to hold all the data compered for further (figure) analysis
        columns = pd.MultiIndex.from_product([major_elements, minor_elements], names=[major, minor])
        data_df = pd.DataFrame(index=df.columns, columns=columns)
        figure = None
        axes = None

        # remove control from the list of elements
        if test not in tests_without_control:
            if len(major_elements) < 2:  # todo: write why it is necessary
                raise Exception('you can not perform this statistical test on one group')
            major_elements.remove(major_control)
        # needs to be after data_df creation so the data frame can have the control as a column

        # data frame to hold the statistical tests results
        stats_df_list = []
        template_stats_df = pd.DataFrame(index=df.columns, columns=['s', 'p', 'p_FDR'])

        # create excel writer to fill up
        excel_path = os.path.join(self.dirs.excels,
                                  'stats {}{} {} {}.xlsx'.format(obj.type, delta_str, test, between))
        excel_writer = pd.ExcelWriter(excel_path)

        # for any major-minor combination
        for major_e in major_elements:
            for minor_e in minor_elements:

                # add a data frame for this specific combination
                stats_df_list.append(template_stats_df.copy())
                if not internal_use:
                    if test not in tests_without_control:
                        stats_df_list[-1].name = '{} {} vs. {} of {} {}'.format(major, major_control, major_e,
                                                                                minor, minor_e)
                        excel_sheet = '{} {} {}'.format(major_control, major_e, minor_e)\
                            .replace('algorithm', 'alg').replace('mediterranean', 'med').replace('fake_time_point', '')
                        # name changes just to limit excel sheet name to the 31 character count limit
                    else:
                        stats_df_list[-1].name = '{} {} of {} {}'.format(major, major_e, minor, minor_e)
                        excel_sheet = '{} {}'.format(major_e, minor_e)
                else:
                    stats_df_list[-1].name = '{}__{}'.format(major_e, minor_e)

                # for each col (bacteria/test/measurement/...)
                for col in df.columns:

                    # define control data if is relevant
                    if test not in tests_without_control:
                        control_data = df[col] \
                            .xs(major_control, level=major, drop_level=False) \
                            .xs(minor_e, level=minor, drop_level=False).dropna()
                        data_df.loc[col, (major_control, minor_e)] = control_data.values
                    else:
                        control_data = None  # just to prevent a warning later on

                    # define the test data no matter the test
                    test_data = df[col] \
                        .xs(major_e, level=major, drop_level=False) \
                        .xs(minor_e, level=minor, drop_level=False).dropna()
                    data_df.loc[col, (major_e, minor_e)] = test_data.values

                    # handling cases where there is missing values in one time_point but not in the other
                    if test in ['wilcoxon', 'ttest_rel']:
                        if between != 'time_point':
                            raise Exception('this test is only relevant when between is time_point')
                        # find the indices that exists in both the control_data and the test_data
                        combined_indices = control_data.index.get_level_values('person').intersection(
                            test_data.index.get_level_values('person'))
                        # filter the data frames to only include these indices
                        control_data = control_data[control_data.index.isin(combined_indices, level='person')]
                        test_data = test_data[test_data.index.isin(combined_indices, level='person')]
                        data_df.loc[col, (major_control, minor_e)] = control_data.values
                        data_df.loc[col, (major_e, minor_e)] = test_data.values

                    # perform the actual test if you have enough samples
                    if test not in tests_without_control and \
                            control_data.unique().shape[0] > minimal_samples and \
                            test_data.unique().shape[0] > minimal_samples:

                        # Compute the Mann-Whitney rank test on samples x and y
                        if test == 'mannwhitneyu':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                mannwhitneyu(control_data, test_data, use_continuity=False, alternative='two-sided')
                            # TODO: think about use_continuity

                        # Calculate the T-test for the means of two independent samples of scores
                        elif test == 'ttest_ind':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                ttest_ind(control_data, test_data, axis=None, equal_var=True, nan_policy='raise')

                        elif test == 'permutation_test':
                            stats_df_list[-1].loc[col, 'p'] = \
                                permutation_test(control_data, test_data, metric=np.average, n_permutations=100)

                        # this is only relevant when between is 'time_point'

                        # Calculates the Wilcoxon signed-rank on TWO RELATED samples of scores, a and b
                        # (non-parametric version of the paired T-test)
                        elif test == 'wilcoxon':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                wilcoxon(control_data, test_data,
                                         zero_method='wilcox', correction=False, alternative='two-sided')
                            # TODO: think about use_continuity (correction) and zero_method

                        # Calculates the T-test on TWO RELATED samples of scores, a and b
                        elif test == 'ttest_rel':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                ttest_rel(control_data, test_data, axis=None)

                        else:
                            raise Exception('test not valid')

                    # tests where test_data is not compared to a control

                    elif test in tests_without_control and \
                            test_data.unique().shape[0] > minimal_samples:

                        if delta is False:
                            raise Exception('this test is only relevant when delta is True')

                        # Calculate the T-test for the mean of ONE group of scores
                        if test == 'ttest_1samp':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                ttest_1samp(test_data, 0, axis=None, nan_policy='raise')

                        # Perform a test that the probability of success is p
                        elif test == 'binom_test':
                            stats_df_list[-1].loc[col, ['s', 'p']] = \
                                binom_test((test_data > 0).sum(), n=test_data.shape[0], p=0.5, alternative='two-sided')
                            # first argument is the number of successes

                    else:
                        # remove in case there are not enough samples
                        stats_df_list[-1].drop(col, inplace=True)
                        # intentionally not done for all columns at once
                        # so problems with test function will be seen as missing values

                # multiple hypothesis correction
                _, stats_df_list[-1]['p_bonferroni'] = bonferroni_correction(stats_df_list[-1]['p'], alpha=self.params.alpha)
                _, stats_df_list[-1]['p_FDR'] = fdr_correction(stats_df_list[-1]['p'], alpha=self.params.alpha)
                stats_df_list[-1].sort_values(f'p_{main_correction_method}', inplace=True)

                # delete from the data_df un-significant results
                significant_indices = stats_df_list[-1]\
                    [stats_df_list[-1][f'p_{main_correction_method}'] < self.params.alpha].index
                non_significant_indices = data_df.index.difference(significant_indices)
                data_df.loc[non_significant_indices, (major_e, minor_e)] = np.nan

                # print for the user
                print('')
                print(test)
                print(stats_df_list[-1].name)
                print('{}/{} {} are significant after bonferroni correction'
                      .format((stats_df_list[-1]['p_bonferroni'] < self.params.alpha).sum(), stats_df_list[-1].shape[0],
                              obj.columns))
                print('{}/{} {} are significant after FDR correction'
                      .format((stats_df_list[-1]['p_FDR'] < self.params.alpha).sum(), stats_df_list[-1].shape[0],
                              obj.columns))
                print('{} {} were not analyzed because they do not have enough samples'
                      .format(df.shape[1] - stats_df_list[-1].shape[0], obj.columns))

                if not internal_use:
                    # excel
                    stats_df_list[-1].to_excel(excel_writer, sheet_name=excel_sheet, freeze_panes=(1, 1))

                    # figure
                    if test not in tests_without_control:
                        figure, axes = \
                            fig_significant_stats(figure_internal=figure, axes_internal=axes,
                                                  curr_data_df=data_df[[major_control, major_e]]
                                                  .xs(minor_e, level=minor, axis=1))
                    else:
                        figure, axes = \
                            fig_significant_stats(figure_internal=figure, axes_internal=axes,
                                                  curr_data_df=data_df
                                                  .xs(major_e, level=major, axis=1, drop_level=False)
                                                  .xs(minor_e, level=minor, axis=1))

        if not internal_use:
            # excel finishes
            excel_writer.save()
            excel_writer.close()
        else:
            return stats_df_list

    def score_models(self, xobj, yobj, cobj, delta=True, minimal_samples=0.9,
                     model_type='linear', n_repeats=1, n_splits=5, random_state=None, hyper_parameters=None,
                     send2queue=False, save=False):
        """
        Calculate n_iterations prediction models or grid search for the hyper parameters for the yobj based on the xobj
        and return the summary statistics on the models scores

        :param xobj: (object) the features object
        :param yobj: (object) the label object
        :param cobj: (object) the covariates object
        :param delta: (bool) whether or not to calculate the models based on the difference between time_point values
        :param minimal_samples: (percentage or int) that needs to be included in the models
        and therefore cannot be deleted in the missing values filtering process
        :param model_type: (str) 'linear' or 'xgb'
        :param n_repeats: (int) number of models to be created for each target
        :param n_splits: (int) number of times to split the data in the KFold or -1 for leave one out
        :param random_state: (int) seed for the random process used for splitting the data
        :param hyper_parameters: (dict) dictionary with hyper parameters for the model {key: [value/s]}
        if value is single will be set directly, if values are multiple will perform a grid search
        :param send2queue: (bool) whether to send the jobs to the queue or run locally
        :param save: (bool) whether or not to save the models

        :return: (pd.DataFrame) scores_df
        """

        # TODO: add model evaluation (ROC, PR, calibration)
        # TODO: add SHAP analysis
        def score_model(x, y, n_splits=1):
            # TODO: not relay on person - maybe x.join(y, how='inner) and then x= joined[:, :-1], y=joined([:, -1])
            # find the indices that exists in the x_df, y_df
            combined_indices = x.dropna().index.get_level_values('person').\
                  intersection(y.dropna().index.get_level_values('person'))

            if len(combined_indices) < minimal_samples:
                print(f'{x_col} {y_col} does not have enough samples')
                return [np.nan]*5

            # filter the data frames to only include these indices
            x = x[x.index.isin(combined_indices, level='person')].sort_values('person')
            y = y[y.index.isin(combined_indices, level='person')].sort_values('person')

            if x.shape[0] != y.shape[0]:
                print(f'{x_col} {y_col} are not the same shape')
                return [np.nan]*5

            # conversion
            x = np.array(x)
            y = np.array(y)

            # check if column is classifiable in order to know which model to use
            if sorted(np.unique(y)) == [0, 1] or sorted(np.unique(y)) == [False, True]:# \
                    # or not all(np.issubdtype(val, np.number) for val in y):
                classifiable = True
            else:
                classifiable = False

            scores = []
            y_pred = []
            y_true = []


            # linear model
            if model_type == 'linear':
                model_instance = LinearRegression(**direct_hyper_params)
                score = 'r2'

            # TODO: add logistic model

            # xgb model
            elif model_type == 'xgb' and classifiable:
                model_instance = XGBClassifier(**direct_hyper_params)
                score = 'roc_auc'
            elif model_type == 'xgb' and not classifiable:
                model_instance = XGBRegressor(objective='reg:squarederror', **direct_hyper_params)
                # objective='reg:squarederror' is the default and is written explicitly just to avoid a warning
                score = 'r2'
            else:
                raise Exception('model not valid')

            # create K folds of the data and do for each fold
            number_of_splits = len(x) if n_splits == -1 else n_splits

            kf = RepeatedKFold(n_splits=number_of_splits, n_repeats=n_repeats, random_state=random_state)  # shuffle=True

            # grid search mode
            if len(search_hyper_params) != 0:

                model = GridSearchCV(model_instance, search_hyper_params, scoring=score, cv=kf)

                x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.5, random_state=random_state)

                model.fit(x_train, y_train)

                if n_splits == -1:
                    y_true.append(y_test)
                    y_pred.append(model.predict(x_test))
                else:
                    scores.append(model.score(x_test, y_test))

            # not grid search mode
            else:

                model = model_instance

                for train_index, test_index in kf.split(x):
                    x_train, x_test = x[train_index], x[test_index]
                    y_train, y_test = y[train_index], y[test_index]

                    model.fit(x_train, y_train)

                    if n_splits == -1:
                        y_true.append(y_test)
                        y_pred.append(model.predict(x_test))
                    else:
                        scores.append(model.score(x_test, y_test))

            if n_splits == -1:
                if score == 'r2':
                    scores = r2_score(np.array(y_true).flatten(), np.array(y_pred).flatten())
                    r_pearson, p_pearson = pearsonr(np.array(y_true).flatten(), np.array(y_pred).flatten())
                    r_spearman, p_spearman = spearmanr(np.array(y_true).flatten(), np.array(y_pred).flatten())
            # save/load the model
            # if save:
            #     pickle.dump(model, open(os.path.join(self.dirs.models,
            #                                          f'model_{xobj.type}_{x_col}_{yobj.type}_{y_col}_' +
            #                                          f'{delta_str}_{model_type}.sav'), 'wb'))
                # model = pickle.load(open('model.sav', 'rb'))

            return [scores, r_pearson, p_pearson, r_spearman, p_spearman]

        # retrieving the data frames from the objects
        if (cobj is None or type(cobj) == self.Object) and type(xobj) == self.Object and type(yobj) == self.Object:
            # TODO: check what happnes if cobj is None
            x_df = cobj.df.join(xobj.df, how='inner') if cobj else xobj.df
            y_df = yobj.df
        else:
            raise Exception('cobj, xobj or yobj type is not object')

        # change the values in the data frames to be the change in value between two time points:
        # any time point and the control time point
        if delta:
            x_df = get_delta_df(x_df, self.params.controls['time_point'])
            y_df = get_delta_df(y_df, self.params.controls['time_point'])
            delta_str = ' delta'  # for file names
        else:
            delta_str = ''  # for file names

        # convert the minimal_samples percentage from the argument to the number of samples
        if minimal_samples < 1:  # meaning minimal_samples is in percentage
            minimal_samples = round(minimal_samples * x_df.shape[0])
            # needs to be after the delta because delta effects the shape

        # split hyper parameters to those given and those that need to be optimized
        if hyper_parameters is None:
            hyper_parameters = {}
        # split the hyper parameters to ones you can directly give the model instance and to ones you need to search for
        direct_hyper_params = {key: value[0] for key, value in hyper_parameters.items() if len(value) == 1}
        search_hyper_params = {key: value for key, value in hyper_parameters.items() if len(value) != 1}

        # random_state is used by the xgb model to split the data the same way in the KFold
        # if random_state is not None:
        #     np.random.seed(random_state)
        # random_state = np.random.randint(100, size=n_repeats)  # the 100 is arbitrary

        # queue
        original_dir = os.getcwd()
        os.chdir(self.dirs.jobs)
        queue = qp if send2queue else fakeqp

        with queue(jobname='model', _delete_csh_withnoerr=True, q=['himem7.q']) as q:

            q.startpermanentrun()
            tkttores = {}

            # create a model for each combination of x an y
            # intentionally taken from object otherwise x includes covariates
            c_cols = cobj.df.columns.tolist() if cobj else []
            for x_col in xobj.df.columns:
                for y_col in yobj.df.columns:
                # send job
                    tkttores[f'{x_col}_{y_col}'] = q.method(score_model, (x_df[c_cols+[x_col]], y_df[[y_col]], n_splits))

            # summarize the scores using the iterations scores distribution if exist
            scores_df = pd.DataFrame(index=tkttores.keys(), columns=['R2', 'r_pearson', 'p_pearson', 'r_spearman', 'p_spearman'])

            for k, v in tkttores.items():
                scores_df.loc[k, ['R2', 'r_pearson', 'p_pearson', 'r_spearman', 'p_spearman']] = q.waitforresult(v)

        os.chdir(original_dir)

        # drop all the un-ran columns
        scores_df = scores_df.dropna(axis=0, how='all')

        # sort values by score
        scores_df = scores_df.sort_values('R2', ascending=False)

        # save the data frame as excel
        # create excel writer to fill up
        excel_path = os.path.join(self.dirs.excels,
                                  'models {} {}{} {}.xlsx'.format(xobj.type, yobj.type, delta_str, model_type))
        excel_writer = pd.ExcelWriter(excel_path)
        scores_df.to_excel(excel_writer, freeze_panes=(1, 1))
        excel_writer.save()
        excel_writer.close()

        return scores_df

    def corr_datasets(self, xobj, yobj, delta=False, group_by=None, minimal_samples=0.1):
        """
        Correlates between the xobj and the yobj overlapping categories.
        The resulting excel is filtered while the figures are not

        :param xobj: (object) some object
        :param yobj: (object) another object
        :param delta: (bool) whether or not to correlate based on the difference between time_point values
        :param group_by: whether to mean multiple values of a category
        - only relevant for figures
        :param minimal_samples: (percentage or int) that needs to be included in order to correlate
        - only relevant for the excel

        :return: None
        """

        def fig_corr_indices(x__df, y__df):

            # unite both data frames into one
            x__df = pd.DataFrame(x__df.stack()).rename(columns={0: xobj.type})
            y__df = pd.DataFrame(y__df.stack()).rename(columns={0: yobj.type})
            df = x__df.join(y__df, how='outer')

            # set all levels to be the same
            df = pd.DataFrame(df.stack()).reset_index()  # to include what used to be columns
            df = df.rename(columns={df.columns[-3]: xobj.columns, df.columns[-2]: 'obj', df.columns[-1]: 'value'})

            # for the title
            if group_by is not None:
                group_by_str = ' grouped by {}'.format(group_by)
            else:
                group_by_str = ''

            # for each parameter
            for param in df.columns[:-1]:  # -1 so to not include 'value'
                indices = list(set(df.columns) - set([param, 'value']))
                param_values = df[param].unique()

                # for each hue
                for hue in indices:
                    param_df = df.pivot_table(index=indices, columns=param, values='value').dropna().reset_index()
                    # defined here and not outside the loop in order to be renewed in case group by changed it

                    # len 2 to be able to plot, shape > 0 to have values, param different than hue for hue to make sense
                    if len(param_values) == 2 and param_df.shape[0] > 0 and param != hue:

                        # group by
                        if group_by is not None:
                            size = param_df.groupby(list(set([group_by, hue])))[param_values]\
                                .apply(lambda group: len(group.dropna())).tolist()
                            param_df = param_df.groupby(list(set([group_by, hue])))[param_values] \
                                .apply(lambda group: group.dropna().mean()).reset_index()
                            # list(set()) is necessary in case group_by == hue and then it cannot be grouped "twice"
                        else:
                            size = None

                        # palette
                        # if all hue values are in colors
                        if len(set(param_df[hue].unique()) - set(self.params.colors.keys())) == 0:
                            palette = self.params.colors
                        else:
                            palette = sns.color_palette("hls", len(param_df[hue].unique()))

                        # plot
                        fig, ax = plt.subplots()

                        sns.scatterplot(x=param_values[0], y=param_values[1], hue=hue, data=param_df,
                                        alpha=0.3, size=size, sizes=(20, 200), palette=palette, ax=ax)

                        # color x and y labels
                        ax.set_xlabel(param_values[0], color=self.params.colors[param_values[0]])
                        ax.set_ylabel(param_values[1], color=self.params.colors[param_values[1]])

                        # add spearman correlation to legend
                        handles, labels = ax.get_legend_handles_labels()
                        for i in np.arange(len(param_df[hue].unique())+1)[1:]:
                            # [1:] so to skip the title but still take all the hues - without the size
                            hue_df = param_df[param_df[hue] == labels[i]]
                            r, p = spearmanr(hue_df[param_values[0]], hue_df[param_values[1]])#, nan_policy='omit')
                            labels[i] = '{}\nr={:.2f}, p={:.2f}'.format(labels[i], r, p)
                        ax.legend(handles, labels)

                        title = 'correlation {}{}{} colored by {}'.format(param, delta_str, group_by_str, hue)
                        plt.suptitle(title)
                        plt.title('each dot represents a {}'.format(list(param_df.columns[:-2]) + [param]))
                        plt.savefig(os.path.join(self.dirs.figs, title))

        # retrieving the data frames from the objects
        if type(xobj) == self.Object and type(yobj) == self.Object:
            x_df = xobj.df
            y_df = yobj.df
        else:
            raise Exception('xobj and yobj type is not object')

        # change the values in the data frames to be the change in value between two time points:
        # any time point and the control time point
        if delta:
            x_df = get_delta_df(x_df, self.params.controls['time_point'])
            y_df = get_delta_df(y_df, self.params.controls['time_point'])
            delta_str = ' delta'  # for file names
        else:
            delta_str = ''  # for file names

        # find the indices and columns that exists in both the x_df and the y_df
        combined_indices = x_df.index.get_level_values('person').intersection(
                           y_df.index.get_level_values('person'))
        combined_columns = x_df.columns & y_df.columns
        # filter the data frames to only include these indices and columns
        x_df = x_df[x_df.index.isin(combined_indices, level='person')][combined_columns]
        y_df = y_df[y_df.index.isin(combined_indices, level='person')][combined_columns]

        # synchronize indices levels (some time one data frame has different indices levels than the other)
        indices_names = list(set(x_df.index.names) & set(y_df.index.names) - set(['sample', 'time']))
        x_indices2remove = list(set(x_df.index.names) - set(indices_names))
        y_indices2remove = list(set(y_df.index.names) - set(indices_names))
        x_df = x_df.droplevel(x_indices2remove)
        y_df = y_df.droplevel(y_indices2remove)

        # convert the minimal_samples percentage from the argument to the number of samples
        if minimal_samples < 1:  # meaning minimal_samples is in percentage
            minimal_samples = round(minimal_samples * len(combined_indices))
            # needs to be after the delta because delta effects the shape

        # plot
        fig_corr_indices(x_df, y_df)

        # TODO: think if the excel portion is still relevant
        # create a data frame to fill with results
        corr_df = pd.DataFrame(index=combined_columns, columns=['rho', 'p', 'p_FDR'])

        # correlate the data sets for each column
        for col in combined_columns:

            # match the x and y values
            data_df = pd.DataFrame(columns=['x', 'y'])
            data_df['x'] = x_df[col]
            data_df['y'] = y_df[col]

            # remove cases where one or both of the datasets are empty
            data_df = data_df.dropna(how='any')

            # make sure you have enough samples and then correlate x and y
            if data_df.shape[0] > minimal_samples:
                corr_df.loc[col, ['rho', 'p']] = spearmanr(data_df['x'], data_df['y'])
                # The p-value roughly indicates the probability of an uncorrelated system producing datasets that have a
                # Spearman correlation at least as extreme as the one computed from these datasets.
                # The p-values are not entirely reliable but are probably reasonable for datasets larger than 500 or so.

        # dropping all the un-ran "columns"
        corr_df = corr_df.dropna(how='all')

        # fdr correction
        _, corr_df['p_FDR'] = fdr_correction(corr_df['p'], alpha=self.params.alpha)
        corr_df.sort_values('p_FDR', inplace=True)

        # print for the user
        print('')
        print('{} and {}{}'.format(xobj.type, yobj.type, delta_str))
        print('{}/{} {} are significant after FDR correction'
              .format((corr_df['p_FDR'] < self.params.alpha).sum(), len(combined_columns), xobj.columns))
        print('{} {} were not analyzed because they do not have enough samples'
              .format(len(combined_columns) - corr_df.shape[0], xobj.columns))

        # save the data frame as excel
        # create excel writer to fill up
        excel_path = os.path.join(self.dirs.excels,
                                  'correlations {} {}{}.xlsx'.format(xobj.type, yobj.type, delta_str))
        excel_writer = pd.ExcelWriter(excel_path)
        corr_df.to_excel(excel_writer, freeze_panes=(1, 1))
        excel_writer.save()
        excel_writer.close()

    def dim_reduction(self, obj, n_pca_comp=5, n_tsne_comp=2, perplexity='auto', pcs2remove=None):
        """
        Reduces the dimensionality of a samples X features data frame
        and save the best 2 components as figures, colored by each index

        :param obj: (Object) with data frame containing samples X features,
        index levels are used to color and create different figures
        :param n_pca_comp: (int) number of desired components from pca
        :param n_tsne_comp: (int) number of desired components from tsne
        :param perplexity: (int) related to the number of nearest neighbors or 'auto'
        :param pcs2remove: (lst) of pcs index to remove from data

        :return: pca_result, tsne_result, df (after removal of specified PCs)
        """

        def fig_best_components():

            # create a separate figure for each index
            for index_name in df.index.names:

                fig, axes = plt.subplots(1, 2, figsize=[rcParams['figure.figsize'][0] * 2,
                                                        rcParams['figure.figsize'][1]])

                # coloring
                hue = df.index.get_level_values(index_name)
                if len(set(hue.unique()) - set(self.params.colors.keys())) == 0:
                    palette = self.params.colors
                else:
                    palette = sns.color_palette("hls", len(hue.unique()))

                # pca
                sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1],
                                hue=hue, palette=palette, s=100, alpha=0.3, ax=axes[0])
                axes[0].set_title('pca')
                axes[0].set_xlabel('pca component 1 - {:.2%}'.format(pca.explained_variance_ratio_[0]))
                axes[0].set_ylabel('pca component 2 - {:.2%}'.format(pca.explained_variance_ratio_[1]))

                # tsne
                sns.scatterplot(x=tsne_result[:, 0], y=tsne_result[:, 1],
                                hue=hue, palette=palette, s=100, alpha=0.3, ax=axes[1])
                axes[1].set_title('tsne')
                axes[1].set_xlabel('tsne component 1')
                axes[1].set_ylabel('tsne component 2')

                title = 'dim reduction {} colored by {}'.format(obj.type, index_name)
                fig.suptitle(title)
                plt.savefig(os.path.join(self.dirs.figs, title))
                plt.close()

                x = df.index.get_level_values(index_name)

                if n_pca_comp <= 5 and x.unique().shape[0] <= 5:

                    fig, axes = plt.subplots(nrows=1, ncols=n_pca_comp,
                                             figsize=[rcParams['figure.figsize'][0] * n_pca_comp,
                                                      rcParams['figure.figsize'][1]])
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.8, wspace=None, hspace=None)

                    for pc in np.arange(n_pca_comp):

                        print(index_name, pc, pca.explained_variance_ratio_[pc])

                        y = pca_result[:, pc]

                        ax = axes[pc]
                        sns.boxplot(x=x, y=y, ax=ax, palette=palette)
                        add_stat_annotation(ax, x=x, y=y, box_pairs=list(combinations(x.unique(), 2)),
                                            test='Mann-Whitney', text_format='star')

                        if x.unique().shape[0] == 2:
                            categories = dict(zip(x.unique().tolist(), [0, 1]))

                            # Correlation
                            r_p, p_p = pearsonr(x.map(categories), y)
                            r_s, p_s = spearmanr(x.map(categories), y)
                            corr_text = 'Correlation: r_p={:.2f}, p_p={:.2f}, r_s={:.2f}, p_s={:.2f}'.\
                                format(r_p, p_p, r_s, p_s)

                            # Logistic Regression
                            model = Logit(endog=x.map(categories).values, exog=y).fit()
                            logit_text = 'Logistic Regression: coef={:.2f}, pval={:.2f}'.\
                                format(model.params[0], model.pvalues[0])
                            text = '\n' + logit_text + '\n' + corr_text

                        else:
                            text = ''

                        ax.set_title('pca component {}{}'.format(pc, text))
                        ax.set_ylabel('pca component {} - {:.2%}'.format(pc, pca.explained_variance_ratio_[pc]))

                        print(text)

                    title = 'PCs {} stratified by {}'.format(obj.type, index_name)
                    fig.suptitle(title)
                    plt.savefig(os.path.join(self.dirs.figs, title))
                    plt.close()

        # retrieving the data frames from the object
        if type(obj) == self.Object:
            df = obj.df
        else:
            raise Exception('obj is not object')

        # pca
        pca = PCA(n_components=n_pca_comp)
        pca_result = pca.fit_transform(df.values)
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        loadings = pd.DataFrame(loadings, index=df.columns)

        # tsne
        if perplexity == 'auto':
            print("Optimizing perplexity")
            for p in range(1, df.shape[0]):#np.arange(1, np.ceil(df.shape[0]**0.5)+1):
                tsne = TSNE(n_components=n_tsne_comp, perplexity=p,
                            method='exact', init='random' if df.shape[0] <= 50 else 'pca')  # needs to be identicalto below
                _ = tsne.fit_transform(df.values)
                if p == 1 or tsne.kl_divergence_ < min_value:
                    perplexity = p
                    min_value = tsne.kl_divergence_
                    print("perplexity ", p, "kl divergence ", tsne.kl_divergence_)
            print("Chose perplexity ", perplexity)

        tsne = TSNE(n_components=n_tsne_comp, perplexity=perplexity,
                    method='exact', init='random' if df.shape[0] <= 50 else 'pca')  # needs to be indentical to above
        tsne_result = tsne.fit_transform(df.values)

        fig_best_components()

        pca_result_in_original_values = pca.inverse_transform(pca_result)
        df = df.subtract(pca_result_in_original_values[:, pcs2remove].sum(axis=1), axis=0)

        return (pca_result, loadings), tsne_result, df

    def fig_snp_heatmap(self, obj, maximal_filling=0.25, minimal_samples=20, species=None, cmap=None, log_colors=False,
                        annotations=None, annotations_cmaps=None, add_hist=True, add_pairs=True, strain_var=None):
        """
        Plot a heatmap based on the SNP dissimilarity data frame for each species

        :param obj: (Object) with SNP dissimilarity data frame
        :param maximal_filling: (float) maximal percentage of samples with missing values to fill with median value
        JUST for clustering, it will not appear in the heatmap!!
        :param minimal_samples: (int) minimal number of samples to have in order to draw
        :param species: (list) species to plot
        :param cmap: (str) dissimilarity color map name
        :param log_colors: (bool) whether to set the heatmap colors to log scale
        :param annotations: (list of str) columns to annotate by
        :param annotations_cmaps: (dict) key; annotation column name, value: color map name
        :param add_hist: (bool) whether to add dissimilarity histograms with colored peaks
        :param add_hist: (bool) whether to add arrows between pairs of samples of the same individual
        :param strain_var: (str) path to strain variability data frame

        :return: None
        """

        full_df = obj.df.copy()
        colors_df = None
        annotations_df = None
        annotations_labels = []
        strain_var = pd.read_hdf(strain_var) if strain_var is not None else None
        pairs = None

        if annotations is not None:
            # create all the annotation data for the colored bars
            annotations_df = full_df
            annotations_df = annotations_df.reset_index()[[s for s in annotations_df.index.names if '1' in s]] \
                .drop_duplicates().set_index('SampleName1')  # take only sample1 metadata
            annotations_df.columns = annotations_df.columns.str.rstrip('1')  # remove 1 suffix from the names
            annotations_df = annotations_df[annotations]  # limit the annotation to these categories
            annotations_df = annotations_df.iloc[~annotations_df.index.duplicated()]  # drop duplicated samples
            annotations_labels = []
            colors_df = annotations_df.copy()
            # replace values with corresponding colors
            for col in colors_df.columns:
                if len(set(colors_df[col]) - set(self.params.colors.keys())) > 0:  # most probably continuous
                    m = cm.ScalarMappable(
                        norm=colors.Normalize(vmin=-colors_df[col].abs().max(), vmax=colors_df[col].abs().max()),
                        cmap=annotations_cmaps[col] if (annotations_cmaps is not None) & (col in annotations_cmaps.keys())
                             else None)
                    colors_df[col] = colors_df[col].apply(m.to_rgba).to_frame()
                else:  # discrete
                    annotations_labels = annotations_labels + list(np.unique(colors_df[col].values))
                    colors_df[col] = colors_df[col].replace(self.params.colors)

        if add_pairs and 'person1' in list(full_df.index.names) \
                     and 'time_point1' in list(full_df.index.names) \
                     and 'group1' in list(full_df.index.names):
            pairs = full_df[(full_df.index.get_level_values('person1') == full_df.index.get_level_values('person2')) &
                            (full_df.index.get_level_values('time_point1') == self.params.controls['time_point']) &
                            (full_df.index.get_level_values('time_point2') != self.params.controls['time_point'])]
            pairs = pairs.reset_index().set_index('Species')[['SampleName1', 'SampleName2', 'group1']]

        if species:
            full_df = full_df[full_df.index.get_level_values('Species').isin(species)]

        # for each species
        for species in full_df.index.get_level_values('Species').unique():

            # create a square data matrix (drop metadata)
            df = full_df.loc[species].reset_index().pivot_table(index='SampleName2',
                                                               columns='SampleName1',
                                                               values='dissimilarity')

            # limit to samples that have at least maximal_filling percentage of existing dissimilarity measurements
            samples_mask = df.columns[df.isna().sum() < df.shape[0] * maximal_filling].values
            df = df.loc[samples_mask, samples_mask]

            if df.shape[0] >= minimal_samples:

                # find the dissimilarities that are still missing
                # and just for the sake of clustering fill them with medians
                na_mask = df.isna()
                df = df.apply(lambda row: row.fillna(row.median()))

                # plotting
                if df.shape[0] > 1:  # because otherwise you cannot cluster

                    df_linkage = linkage(squareform(df, force='tovector', checks=False),
                                         method='average', metric='euclidean', optimal_ordering=True)
                    # TODO: consider method='max' so to be like in assemblies

                    # necessary to define explicitly otherwise there are differences between the rows and columns

                    # A (n-1) by 4 matrix Z is returned. At the i-th iteration,
                    # clusters with indices Z[i, 0] and Z[i, 1] are combined to form cluster n + i.
                    # A cluster with an index less than n corresponds to one of the original observations.
                    # The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2].
                    # The fourth value Z[i, 3] represents the number of original observations in the newly
                    # formed cluster.

                    # the clustermap and the heatmap are separated because the clustermap does not pass the norm
                    # to the heatmap well as it should
                    # https://github.com/mwaskom/seaborn/pull/1830/files

                    if strain_var is not None:
                        m = cm.ScalarMappable(
                            norm=mpl.colors.Normalize(vmin=0, vmax=strain_var[species].max()), cmap=cm.binary)
                        species_colors_df = strain_var[species].map(m.to_rgba).to_frame().rename(columns={species: ' '})
                        if colors_df is not None:
                            species_colors_df = colors_df.join(species_colors_df)
                    else:
                        species_colors_df = colors_df

                    g = sns.clustermap(df, mask=(df >= 0),  # the mask here is for everything
                                       xticklabels=False, yticklabels=False,
                                       row_linkage=df_linkage, col_linkage=df_linkage,
                                       row_colors=species_colors_df, col_colors=species_colors_df)

                    norm = LogNorm(g.data2d.min().min(), g.data2d.max().max(), clip=False) if log_colors else None

                    sns.heatmap(g.data2d, mask=na_mask.loc[g.data2d.index, g.data2d.columns],
                                xticklabels=False, yticklabels=False,
                                ax=g.ax_heatmap, cbar_ax=g.cax,
                                cmap=cmap, norm=norm,
                                cbar_kws={'label': 'dissimilarity', 'orientation': 'horizontal', 'ticks': None})
                    # cbar ticks are None because otherwise in log scale there are no ticks

                    if pairs is not None:
                        pairs2plot = pairs.loc[species]
                        pairs2plot = pairs2plot[(pairs2plot['SampleName1'].isin(samples_mask)) &
                                                (pairs2plot['SampleName2'].isin(samples_mask))]

                        for _, (s1, s2, group) in pairs2plot.iterrows():
                            s1i = g.data2d.index.tolist().index(s1)
                            s2i = g.data2d.index.tolist().index(s2)
                            g.ax_heatmap.arrow(x=s1i, y=s1i, dx=s2i-s1i, dy=0,
                                               color=self.params.colors[group], head_width=3, head_length=2)

                    expanded_labels = []
                    is_significant = []

                    # statistical test
                    if annotations:
                        stats_df = annotations_df.loc[df.iloc[g.dendrogram_col.reordered_ind].index]
                        stats_df['rank'] = np.arange(stats_df.shape[0])

                        for anno in annotations:

                            # statistics
                            if len(np.unique(stats_df[anno])) == 2:
                                r, p = mannwhitneyu(
                                    x=stats_df.loc[stats_df[anno] == np.unique(stats_df[anno])[0], 'rank'].values,
                                    y=stats_df.loc[stats_df[anno] == np.unique(stats_df[anno])[1], 'rank'].values,
                                    use_continuity=True, alternative='two-sided')  # TODO: think about use_continuity
                            else:
                                r, p = spearmanr(stats_df[anno], stats_df['rank'])
                            r = round(r, 2)
                            p = round(p, 2)

                            # rand index between the annotation and the clusters
                            org_clusters = annotations_df.loc[df.index, anno].values.flatten()
                            new_clusters = fcluster(df_linkage, t=len(np.unique(org_clusters)),
                                                    criterion='maxclust')  # TODO: optimize parameters

                            RI = round(adjusted_rand_score(org_clusters, new_clusters), 2)  # TODO: optimize parameters

                            # TODO: think what should be the t for inconsistent fcluster
                            # new_clusters2 = fcluster(df_linkage, t=1, criterion='inconsistent')
                            # rand_index2 = adjusted_rand_score(org_clusters, new_clusters2)

                            expanded_labels.append('{} (r={}, p={}, RI={})'.format(anno, r, p, RI))
                            is_significant.append(p)

                    if strain_var is not None:
                        expanded_labels.append(f'variable positions (max {int(strain_var[species].max() * 100)}%)')
                    if len(expanded_labels) > 0:
                        g.ax_col_colors.set_yticklabels(expanded_labels)

                    # title
                    title = '{}\n{}\n{}'.format(obj.type, segata_name(species), species)
                    g.fig.suptitle(title)

                    # annotations legend
                    for label in annotations_labels:
                        g.ax_col_dendrogram.bar(0, 0, color=self.params.colors[label], label=label, linewidth=0)

                    # position of plot
                    plt.subplots_adjust(left=0.02, right=0.76, bottom=0.11, top=0.85)

                    # position of annotations legend
                    g.ax_col_dendrogram.legend(ncol=1, frameon=False, bbox_to_anchor=(-0.035, 1))

                    # position of dissimilarity legend
                    g.cax.set_position([.23, .06, .5, .01])

                    # dissimilarity histograms with colored peaks
                    if add_hist:
                        diss_tril = g.data2d[~na_mask.loc[g.data2d.index, g.data2d.columns]].values
                        diss_tril = diss_tril[np.tril_indices(n=diss_tril.shape[0], k=-1)]

                        special_ax = inset_axes(g.ax_heatmap, width="40%", height=1.0, loc='lower left')
                        bins_freq, _, patches = special_ax.hist(diss_tril, bins='auto')
                        peaks_bin, _ = find_peaks(bins_freq)  # TODO: optimize
                        [patches[bin].set_color('orange') for bin in peaks_bin]
                        special_ax.set_title('dissimilarity distribution', color='white')
                        special_ax.set_xticklabels([])
                        special_ax.set_yticklabels([])

                    if not os.path.exists(os.path.join(self.dirs.figs, obj.type)):
                        os.makedirs(os.path.join(self.dirs.figs, obj.type))
                    plt.savefig(os.path.join(self.dirs.figs, obj.type, species), pad_inches=0.5)

                    if any([p < self.params.alpha for p in is_significant]):
                        if not os.path.exists(os.path.join(self.dirs.figs, obj.type, 'significant')):
                            os.makedirs(os.path.join(self.dirs.figs, obj.type, 'significant'))
                        plt.savefig(os.path.join(self.dirs.figs, obj.type, 'significant', species), pad_inches=0.5)

                    plt.close()

    def fig_snp_scatter_box(self, obj, subplot='group', subplot_order=None,
                            minimal_between_comparisons=45, minimal_within_comparisons=10,
                            species=25, height=12, aspect=0.5, whis=[5, 95],
                            sum_not_mean=True):
        """
        Plot a the dissimilarity distribution between and within people based on the SNP dissimilarity data frame
        for each species

        :param obj: (Object) with SNP dissimilarity data frame
        :param subplot: (str) group or time_point as subplots
        :param minimal_between_comparisons: (int) minimal number of comparisons between people per specie
        :param minimal_within_comparisons: (int) minimal number of comparisons within people per specie
        :param species: (int or list) int of number of species per plot or a list of specific species to plot
        :param height: (int) sns.FacetGrid height argument
        :param aspect: (float) sns.FacetGrid aspect argument
        :param whis: () sns.bosxplot whis argument
        :param sum_not_mean: (bool) sum (True) or mean (False) replacements

        :return: None
        """

        def scatter_box_plot(**kwargs):
            # retract the data
            sub_data = kwargs.pop('data')

            if not summary:

                # between individuals
                if not sub_data.loc[~same_person & same_minor].empty:
                    sns.boxplot('dissimilarity', 'Species', minor,
                                data=sub_data.loc[~same_person & same_minor],
                                hue_order=np.unique(sub_data.loc[~same_person & same_minor, minor]),
                                width=0.5, whis=whis, showfliers=False, **kwargs)

                # within individuals
                if not sub_data.loc[same_person].empty:
                    hue = major if not all(sub_data.loc[same_person, 'person1'].isin(list(self.params.colors.keys()))) else 'person1'
                    p = sns.scatterplot('dissimilarity', 'Species', hue,
                                    data=sub_data.loc[same_person],
                                    hue_order=np.unique(sub_data.loc[same_person, hue]),
                                    alpha=0.2, legend='brief', **kwargs)

                    # replacements
                    if 'replacement' in sub_data.columns:
                        if sum_not_mean:
                            n_replacements = sub_data.loc[same_person].groupby(['Species', 'person1'])['replacement'].any()\
                                .groupby('Species').sum().replace(False, 0).replace(True, 1)
                        else:
                            n_replacements = sub_data.loc[same_person].groupby(['Species', 'person1'])['replacement'].any() \
                                .replace(False, 0).replace(True, 1).groupby('Species').mean()
                        for s in n_replacements.index:
                            p.text(1/20, s, f'{n_replacements.loc[s]:.0%}')

                # between individuals - significance
                if minimal_between_comparisons is not None:
                    sns.scatterplot(1/15, 'Species', minor,
                                    data=stats_between_df.loc[stats_between_df[major].isin(np.unique(sub_data[major])) &
                                                              stats_between_df['Species'].isin(np.unique(sub_data['Species']))],
                                    hue_order=np.unique(stats_between_df.loc[
                                                            stats_between_df[major].isin(np.unique(sub_data[major])) &
                                                            stats_between_df['Species'].isin(np.unique(sub_data['Species'])),
                                                            minor]),
                                    marker='*', s=120, legend=False, **kwargs)

                # within individuals - significance
                if minimal_within_comparisons is not None:
                    sns.scatterplot(1/20, 'Species', major,
                                    data=stats_within_df.loc[stats_within_df[major].isin(np.unique(sub_data[major])) &
                                                             stats_within_df['Species'].isin(np.unique(sub_data['Species']))],
                                    hue_order=np.unique(stats_within_df.loc[
                                                            stats_within_df[major].isin(np.unique(sub_data[major])) &
                                                            stats_within_df['Species'].isin(np.unique(sub_data['Species'])),
                                                            major]),
                                    marker='*', s=120, legend=False, **kwargs)

            else:

                # within individuals
                if not sub_data.loc[same_person].empty:
                    matrix = sub_data.loc[same_person].copy()
                    matrix = matrix.groupby(['Species', 'person1'])['replacement'].any().unstack().fillna(False)
                    matrix = matrix.apply(lambda col: col.replace(True, pd.to_numeric(col.name, errors='ignore')),
                                          axis=0).replace(False, 0)
                    matrix.loc[:, set(people_order) - set(matrix.columns)] = 0
                    matrix = matrix.loc[:, people_order]

                    # TODO: bug and a problem if the person1 identifiers are not sequential and in equal distance
                    cmap = None if not all(sub_data.loc[same_person, 'person1'].isin(list(self.params.colors.keys()))) \
                        else ['white'] + [self.params.colors[str(k)] for k in matrix.columns]
                    h = sns.heatmap(matrix, linewidths=1, linecolor='lightgrey', cbar=False, cmap=cmap)

                    idx = np.unique(sub_data.loc[same_person, major])[0]  # assuming there aren't any bugs with major
                    h.set_xticklabels([f'{tick.get_text()}\n{people_info.loc[(tick.get_text(), idx)]}'
                                       if (tick.get_text(), idx) in people_info.index else f'{tick.get_text()}'
                                       for tick in h.get_xticklabels()])

        # data manipulation
        df = obj.df
        df = df.reset_index()

        # species argument
        if type(species) == list:  # by specific species
            df = df[df['Species'].isin(species)]
            species = len(species)

        df = df[df['SampleName1'] != df['SampleName2']]
        df['s1'] = np.minimum(df['SampleName1'], df['SampleName2'])
        df['s2'] = np.maximum(df['SampleName1'], df['SampleName2'])
        df = df.drop_duplicates(['Species', 's1', 's2'])

        # can subplot groups (control/test1/test2/...) or time points (0/1/2...)
        subplot_options = ['group', 'time_point']  # this list has to be length 2
        if subplot not in subplot_options:
            raise Exception('invalid between')

        # names of elements
        major = subplot
        minor = subplot_options[abs(subplot_options.index(major) - 1)]  # the index of the not major element

        # it does not need to be redefined after filtering the data frame because it uses the indices
        same_person = df['person1'] == df['person2']
        same_major = df['{}1'.format(major)] == df['{}2'.format(major)]
        same_minor = df['{}1'.format(minor)] == df['{}2'.format(minor)]

        df.loc[~same_major, major] = \
            np.minimum(df['{}1'.format(major)], df['{}2'.format(major)]) + ' vs. ' + \
            np.maximum(df['{}1'.format(major)], df['{}2'.format(major)])
        df.loc[same_major, major] = df['{}1'.format(major)]

        df.loc[~same_minor, minor] = \
            np.minimum(df['{}1'.format(minor)], df['{}2'.format(minor)]) + ' vs. ' + \
            np.maximum(df['{}1'.format(minor)], df['{}2'.format(minor)])
        df.loc[same_minor, minor] = df['{}1'.format(minor)]

        # statistics
        # TODO: fix bug in saving excel
        # between individuals
        if minimal_between_comparisons is not None:
            obj4stats_between = Study.Object(obj_type='{}_'.format(obj.type), columns='bacteria')
            obj4stats_between.df = df.loc[~same_person & same_minor].set_index(['Species', major, minor], append=True)
            obj4stats_between.df = obj4stats_between.df['dissimilarity'].unstack('Species')

            if len(obj4stats_between.df.index.get_level_values(minor).unique()) > 1:
                # this is a test to check if we randomize the labels we still get so many significant results
                # obj4stats_between.df = obj4stats_between.df.reset_index(minor)
                # for curr_major in np.unique(obj4stats_between.df.index.get_level_values(major)):
                #     obj4stats_between.df.loc[obj4stats_between.df.index.get_level_values(major) == curr_major, minor] = \
                #         shuffle(obj4stats_between.df.loc[obj4stats_between.df.index.get_level_values(major) == curr_major, minor].values)  # from sklearn.utils import shuffle
                # obj4stats_between.df = obj4stats_between.df.set_index(minor, append=True)
                stats_between_df = self.comp_stats(obj4stats_between, test='mannwhitneyu', between=minor,
                                                   minimal_samples=minimal_between_comparisons, internal_use=True)

                for i in np.arange(len(stats_between_df)):
                    stats_between_df[i][major] = stats_between_df[i].name.split('__')[1]
                    stats_between_df[i][minor] = stats_between_df[i].name.split('__')[0]
                stats_between_df = pd.concat(stats_between_df)

                stats_between_df = stats_between_df[stats_between_df['p_FDR'] < self.params.alpha].reset_index()

            else:
                stats_between_df = pd.DataFrame(columns=['Species', major, minor])

        # within individuals
        if minimal_within_comparisons is not None:
            obj4stats_within = Study.Object(obj_type='{}_'.format(obj.type), columns='bacteria')
            obj4stats_within.df = df.loc[same_person].set_index(['Species', major, minor], append=True)
            obj4stats_within.df = obj4stats_within.df['dissimilarity'].unstack('Species')

            if len(obj4stats_within.df.index.get_level_values(major).unique()) > 1:
                stats_within_df = self.comp_stats(obj4stats_within, test='mannwhitneyu', between=major,
                                                  minimal_samples=minimal_within_comparisons, internal_use=True)

                for i in np.arange(len(stats_within_df)):
                    stats_within_df[i][major] = stats_within_df[i].name.split('__')[0]
                    stats_within_df[i][minor] = stats_within_df[i].name.split('__')[1]
                stats_within_df = pd.concat(stats_within_df)

                stats_within_df = stats_within_df[stats_within_df['p_FDR'] < self.params.alpha].reset_index()

            else:
                stats_within_df = pd.DataFrame(columns=['Species', major, minor])

        # removing species that lack inter or intra person data
        df = df[df['Species'].isin(list(
            set(df.loc[same_person, 'Species']).intersection(set(df.loc[~same_person & same_minor, 'Species']))))]

        if minimal_between_comparisons is not None and minimal_within_comparisons is not None:
            combined_species = set(stats_between_df['Species']).union(set(stats_within_df['Species']))  # significant
            combined_species = combined_species.intersection(set(df['Species']))
            df = df[df['Species'].isin(combined_species)]
            stats_between_df = stats_between_df[stats_between_df['Species'].isin(combined_species)]
            stats_within_df = stats_within_df[stats_within_df['Species'].isin(combined_species)]

        # sorts species and people df by the number of replacements
        if 'replacement' in df.columns:
            df['replacement'] = df['replacement'].map({True: 1, False: 0})
            if sum_not_mean:
                species_order = list(df.loc[same_person].groupby(['Species', 'person1'])['replacement'].any().
                                     groupby('Species').sum().sort_values(ascending=False).index)
            else:
                species_order = list(df.loc[same_person].groupby(['Species', 'person1'])['replacement'].any() \
                                .replace(False, 0).replace(True, 1).groupby('Species').mean().sort_values(ascending=False).index)

        else:
            species_order = list(df.loc[same_person].groupby('Species').apply(len).sort_values(ascending=False).index)

        # this is just to know who are all the possible people in the same person comparisons
        str_values = []
        int_values = []
        for v in df.loc[same_person, 'person1'].unique():
            if v.isnumeric():
                int_values.append(int(v))
            else:
                str_values.append(v)
        people_order = sorted(str_values) + [str(v) for v in sorted(int_values)]

        title = '{} distribution by {}'.format(obj.type, subplot)

        people_info = df.loc[same_person].groupby(['Species', 'person1', 'group1'])['replacement'].any()
        people_info.reset_index().to_excel(os.path.join(self.dirs.excels, title + '.xlsx'))
        people_info = people_info.groupby(['person1', 'group1']).apply(lambda g: f'({g.sum()}/{g.shape[0]})')

        df['Species'] = df['Species'].astype('category').cat.set_categories(species_order)
        df = df.sort_values(['Species'])
        df['Species'] = df['Species'].astype(str)
        # for unexplainable reason this effects the replacements decimal point in the figure
        # if not returned to type str will show all species in all figures

        # x axis limits
        max_x = 10**-1
        min_x = self.params.dissimilarity_threshold/3
        df['dissimilarity'] = df['dissimilarity'].clip(lower=self.params.dissimilarity_threshold, upper=max_x)

        # split to multiple figures
        for i_figure in np.arange(0, len(np.unique(df['Species'])), species):

            sub_df = df[df['Species'].isin(species_order[i_figure:i_figure+species])]
            sub_df['Species'] = sub_df['Species'].astype('category').cat.set_categories(species_order[i_figure:i_figure+species])
            # if not type categorial there is a bug in y values order (https://github.com/mwaskom/seaborn/issues/1306)

            for summary in [False, True]:

                # plotting
                g = sns.FacetGrid(sub_df, col=major,
                                  col_order=np.unique(df[major]) if subplot_order is None else subplot_order,
                                  height=height, aspect=aspect, dropna=False, sharex=not summary)
                g = g.map_dataframe(scatter_box_plot, palette=self.params.colors)

                if not summary:
                    plt.xscale('log')
                    plt.xlim([min_x, max_x])

                # ticks
                ax = g.axes.flatten()[0]
                ax.set_yticklabels([f'{segata_name(tick.get_text())} ({tick.get_text()})'
                                    for tick in ax.get_yticklabels()])

                # titles and labels
                g.set_axis_labels(x_var='dissimilarity' if not summary else
                                        'person\n(# replaced species / # compared species)', y_var='species')
                g.set_titles(row_template='{row_name}' if not summary else None, col_template='{col_name}')
                for ax in g.axes.flatten():
                    if ax.get_title() in self.params.colors.keys():
                        color = self.params.colors[ax.get_title()]
                    else:
                        color = 'black'
                    ax.set_title(label=ax.get_title(), color=color)

                plt.suptitle(title)

                # figure
                plt.tight_layout()
                plt.subplots_adjust(top=0.88)

                # saving
                title4saving = '{}{}s{}-{}'.format(title, ' ' if not summary else ' summary ', i_figure, i_figure+species)
                if not os.path.exists(os.path.join(self.dirs.figs, obj.type)):
                    os.makedirs(os.path.join(self.dirs.figs, obj.type))
                plt.savefig(os.path.join(self.dirs.figs, obj.type, title4saving))
                plt.close()

        return species_order

    def fig_strain_replacements(self, obj, sum_not_mean=True, category='person1', hue='group1'):
        """
        Sum/mean and plot (histogram) replacements per category

        :param obj: (Object) with SNP dissimilarity data frame and replacements
        :param sum_not_mean: (bool) sum (True) or mean (False) replacements
        :param category: (str) index level name to count by
        :param hue: (str) index level name to color by

        :return: n_replacements (pd.DataFrame)
        """
        df = obj.df
        df = df.reset_index()
        df = df[df['SampleName1'] != df['SampleName2']]
        df['s1'] = np.minimum(df['SampleName1'], df['SampleName2'])
        df['s2'] = np.maximum(df['SampleName1'], df['SampleName2'])
        df = df.drop_duplicates(['Species', 's1', 's2'])
        df = df.dropna(subset=['replacement'])
        df['replacement'] = df['replacement'].map({True: 1, False: 0})

        groups = [category, hue] if hue is not None else [category]
        if sum_not_mean:
            replacements = df.groupby(groups)['replacement'].sum().sort_values(ascending=False)
        else:
            replacements = df.groupby(groups)['replacement'].mean().sort_values(ascending=False)

        if category[-1] in ['1', '2']:
            category = category[:-1]
        category = category.lower()
        if hue is not None:
            if hue[-1] in ['1', '2']:
                hue = hue[:-1]
            hue = hue.lower()
            replacements.index.names = [category, hue]
        else:
            replacements.index.name = category

        # histogram
        g = sns.FacetGrid(replacements.reset_index(), hue=hue, palette=self.params.colors,
                          sharex=False, sharey=False, margin_titles=True)
        g = g.map(sns.distplot, 'replacement', hist=True, kde=False, norm_hist=True)  # norm_hist=False does not work as expected
        g = g.add_legend()

        obj_type = obj.type.split(' ')[0].lower()
        obj_type = obj_type + ' ' if obj_type in ['oral', 'gut'] else ''
        plt.title('{}strain replacements per {} histogram'.format(obj_type, category))
        if sum_not_mean:
            plt.xlabel('number of replacements')
        else:
            plt.xlabel('mean replacements')
        plt.ylabel('number of {}'.format(category))

        # statistical test
        hues = replacements.index.get_level_values(hue).unique()
        if len(hues) == 2:
            _, p = mannwhitneyu(x=replacements.xs(hues[0], level=hue).dropna().tolist(),
                                y=replacements.xs(hues[1], level=hue).dropna().tolist(),
                                use_continuity=True, alternative="two-sided", axis=0, method="auto")
            plt.text(x=0.9, y=0.9, s=f'p={p:.2f}', transform=g.ax.transAxes)

        plt.savefig(os.path.join(self.dirs.figs, 'replacements per {}{} {}'.format(obj_type, category, 'sum' if sum_not_mean else 'mean')))

        return replacements

    def fig_correlation_heatmap(self, xobj, yobj, minimal_samples=20):
        """
        Correlate all columns between xobj and yobj and plot the results in a heatmap

        :param xobj: (Object) to plot on the x-axis
        :param yobj: (Object) to plot on the y-axis
        :param minimal_samples: (int) to have in order to check correlations

        :return: None
        """

        x_df = xobj.df
        y_df = yobj.df

        # clean abundance species names
        if 'abundance' in xobj.type:
            x_df.columns = [col.split('|')[-1].replace('sSGB_', 'SGB') for col in x_df.columns]

        # manipulate dissimilarity data frame
        if 'dissimilarity' in x_df.columns:
            # to only contain self comparison
            x_df = x_df[x_df.index.get_level_values('person1') == x_df.index.get_level_values('person2')]
            # to match non-dissimilarity data frames
            x_df['person'] = x_df.index.get_level_values('person1')
            x_df = x_df.set_index('person', append=True)
            # transform data frame so Species would be columns
            x_df = x_df.pivot_table(index='person', columns='Species', values='dissimilarity')

        # sync data frames
        combined_indices_levels = set(x_df.index.names) & set(y_df.index.names)
        x_df.index = x_df.index.droplevel(list(set(x_df.index.names) - combined_indices_levels))
        y_df.index = y_df.index.droplevel(list(set(y_df.index.names) - combined_indices_levels))

        combined_indices_values = x_df.index.intersection(y_df.index)
        x_df = x_df.loc[combined_indices_values]
        y_df = y_df.loc[combined_indices_values]

        # calculate correlations
        corr_df = pd.DataFrame(index=x_df.columns, columns=y_df.columns)
        pval_df = pd.DataFrame(index=x_df.columns, columns=y_df.columns)

        for x_col in x_df.columns:
            for y_col in y_df.columns:

                # sync columns
                x = x_df[x_col].dropna()
                y = y_df[y_col].dropna()
                combined_indices = x.index.intersection(y.index)
                x = x.loc[combined_indices]
                y = y.loc[combined_indices]

                if x.shape[0] > minimal_samples:
                    corr_df.loc[x_col, y_col], pval_df.loc[x_col, y_col] = pearsonr(x, y)
                    # TODO: think if should be spearman
        # TODO: multiple comparisons correction

        # excel
        results_df = pd.DataFrame()
        results_df['r'] = corr_df.values.flatten().dropna(how='all')
        results_df['p'] = pval_df.values.flatten().dropna(how='all')
        results_df.to_csv(os.path.join(self.dirs.excels, 'correlation {} {}.csv'
                                       .format(xobj.type.replace('_', ' '), yobj.type.replace('_', ' '))))
        # TODO: us excel writer like everywhere else

        # plot
        pval_df = pval_df.dropna(how='all', axis=0)
        corr_df = corr_df.loc[pval_df.index, pval_df.columns]
        fig, ax = plt.subplots(figsize=(rcParams['figure.figsize'][0]*2, rcParams['figure.figsize'][1]))
        ax = sns.heatmap(corr_df.T.astype(float), mask=(pval_df.T < self.params.alpha), cmap='coolwarm')
        # TODO: replace with clustermap
        plt.title('Correlation between {} and {} {}'.format(xobj.type.replace('_', ' '),
                                                            yobj.type.replace('_', ' '), yobj.columns))
        plt.xlabel(xobj.type.replace('_', ' '))
        plt.ylabel(yobj.type.replace('_', ' '))

        plt.savefig(os.path.join(self.dirs.figs, 'correlation {} {}'
                                 .format(xobj.type.replace('_', ' '), yobj.type.replace('_', ' '))))


# general file handling
def ftp_download(address, username, password, directories, skip_files=None, destination=os.getcwd(), check_only=True):
        """
        Download files from a FTP address not including those in skip_files to the destination directory

        :param address: (str) FTP address to connect to
        :param username: (str) username to use in the FTP server
        :param password: (str) password to use in the FTP server
        :param directories: (str, list) directories to download from the FTP server
        :param skip_files: (str, list) files and sub-directories that should not be downloaded
        :param destination: (str) directory for the downloaded files to be saved in
        :param check_only: (bool) only compare the files in the address against those in the destination
        instead of downloading them

        :return: None
        """

        # connecting
        ftp = FTP(address, username, password)  # using the default port

        for Dir in directories:
            files_in_dir = ftp.nlst(Dir)  # list directory content

            for ftp_file_name in files_in_dir:
                linux_file_name = os.path.join(destination, os.path.basename(ftp_file_name))

                if ftp_file_name not in skip_files:

                    # check file existence and size
                    if check_only:
                        ftp.sendcmd('TYPE I')  # binary mode - necessary for size check
                        # file existence
                        if not os.path.exists(linux_file_name):
                            raise Warning('missing file: {}'.format(linux_file_name))
                        # file size
                        elif ftp.size(ftp_file_name) != os.path.getsize(linux_file_name):
                            raise Warning('different file size: {}'.format(linux_file_name))

                    # download the file
                    elif not os.path.exists(linux_file_name):
                        ftp.retrbinary('RETR ' + ftp_file_name, open(linux_file_name, 'wb').write)

        ftp.quit()


def decompress_files(file_type, input_dir=os.getcwd(), output_dir=os.getcwd(), regex_name='*', delete=False):
        """
        Decompress files of type file_type that matches the regex_name from the input_dir to the output_dir
        and delete the compressed files according to "delete"

        :param file_type: (str) 'tar'
        :param input_dir: (str) input directory to look for files
        :param output_dir: (str) output directory to save decompressed files
        :param regex_name: (str) limit the files to decompress to those that match this name pattern
        (do not include file type)
        :param delete: (bool) if files should be deleted at the end of the process

        :return: None
        """

        compressed_files = glob(input_dir + '/' + regex_name + '.' + file_type)

        for file_name in compressed_files:

            # decompressing
            if file_type == 'tar':
                file = tarfile.open(file_name)
                file.extractall(output_dir)
                file.close()
            else:
                raise Exception('file_type is not implemented in this function')

            # deleting the compressed files
            if delete:
                os.remove(file_name)


# data frame computations
def get_delta_df(regular_df, control_time_point, divide=False):
    """
    Subtract from each time point values the self.params.control_time values

    :param regular_df: (pd.DataFrame) data frame to calculate for the delta
    :param control_time_point: time point values to subtract for all other time points values
    :param divide: (bool) default False mean subtract value True means divide

    :return: delta_df (pd.DataFrame)
    """

    # list of all indices columns TO REMOVE in order to not have a contradicting indices
    # between different time points
    index_names2remove = list(regular_df.index.names)
    index_names2remove.remove('person')
    index_names2remove.remove('time_point')

    # copy of the data frame in order to not delete the additional indices from the returned data frame
    delta_df = regular_df.copy()
    regular_df = regular_df.droplevel(index_names2remove)

    # all possible time points
    time_points = delta_df.index.get_level_values('time_point').unique().to_list()
    time_points.remove(control_time_point)

    # subtract from each time point values the control time values
    for time_point in time_points:
        idx = delta_df.xs(time_point, level='time_point', drop_level=False).index
        idx_person = idx.get_level_values('person')

        curr_values = regular_df.xs(time_point, level='time_point')
        control_values = regular_df.xs(control_time_point, level='time_point')

        if divide:
            delta_df.loc[idx] = (curr_values / control_values).loc[idx_person].values
            sign = '/'
        else:
            delta_df.loc[idx] = (curr_values - control_values).loc[idx_person].values
            sign = '-'

        delta_df = delta_df.rename(
            index={time_point: '{}{}{}'.format(time_point, sign, control_time_point)},
            level='time_point')

    # delete the control time points
    idx = delta_df.xs(control_time_point, level='time_point', drop_level=False).index
    delta_df = delta_df.drop(idx)

    return delta_df


def time_series2time_point(person_df, days_between_time_points=6):
    """
    Takes the time index column and adds to it time_point columns for each group of samples that have less
    than days_between_time_points between them

    :param person_df: (pd.DataFrame) single person data frame to apply the function on
    :param days_between_time_points: (int) maximum amount of days that can be between samples that will still be
    considered as a single time_point

    :return: person_df (pd.DataFrame)
    """

    person_df['diff'] = person_df.index.get_level_values('time')
    person_df = person_df.sort_values('diff')  # which is actually 'time'
    person_df['diff'] = person_df['diff'].diff().apply(lambda x: x.days)

    break_points = person_df['diff'][person_df['diff'] > days_between_time_points].index
    break_points = break_points.union([person_df.index[-1]])  # add last row

    person_df['time_point'] = np.nan
    previous_break_point = person_df.index[0]
    for time_point, current_break_point in enumerate(break_points):
        person_df.loc[previous_break_point:current_break_point, 'time_point'] = time_point
        previous_break_point = current_break_point

    person_df = person_df.drop('diff', axis=1).set_index('time_point', append=True)

    return person_df


def add_cgm(samples_df, cgm_df, delta_days=-7, glucose_threshold=140):
    """
    Adds percentage of cgm measurements above glucose_threshold from these amount of days_back
    to each sample in samples_df

    :param samples_df: (pd.DataFrame) to add to each sample the relevant measurements
    :param cgm_df: (pd.DataFrame) to take the relevant cgm measurements from
    :param delta_days: (int) cgm time in days to consider before (negative) or after (positive) sample time
    :param glucose_threshold: (int) to count measurements above it

    :return: samples_df (pd.DataFrame)
    """

    # TODO: consider the fact that different cgms have different baseline values
    samples_df['time_above{}'.format(glucose_threshold)] = np.nan
    for i, sample in samples_df.reset_index().iterrows():

        if sample['person'] in cgm_df.index.get_level_values('person'):
            curr_cgm = cgm_df.xs(sample['person'], level='person')  # this person
            delta_time = (curr_cgm.index.get_level_values('time').tz_convert('UTC') -
                          sample['time']).days  # conversion to UTC should be part of the LabData

            if delta_days < 0:
                curr_cgm = curr_cgm[(delta_days <= delta_time) & (delta_time <= 0)]  # these days back
            else:
                curr_cgm = curr_cgm[(0 <= delta_time) & (delta_time <= delta_days)]  # these days forward

            if curr_cgm.shape[0] != 0:
                samples_df.iloc[i, -1] = (curr_cgm['GlucoseValue'] > glucose_threshold).sum() / curr_cgm.shape[0]
                # assuming the last column is the column added

    return samples_df


# microbiome abundance
def get_diversity_df(abundance_df):
    """
    Compute the Shanon's alpha diversity index (using log2) based on an "abundance_df"

    :param abundance_df: (pd.DataFrame) data frame to compute on

    :return: (pd.DataFrame) "diversity_df"
    (same as abundance_df just instead of bacteria columns there is a single diversity column)
    """

    # Shannon's alpha diversity index = -sum(Pi*log2(Pi))

    # revert values to their pre log state
    if not (0 <= abundance_df.min().min() and abundance_df.max().max() <= 1):
        abundance_df = (10 ** abundance_df)

    diversity_df = pd.DataFrame(-(abundance_df * np.log2(abundance_df)).sum(axis=1))
    diversity_df.columns = ['diversity']

    return diversity_df


def get_richness_df(abundance_df):
    """
    Compute the richness based on an "abundance_df"

    :param abundance_df: (pd.DataFrame) data frame to compute on

    :return: (pd.DataFrame) "richness_df"
    (same as abundance_df just instead of bacteria columns there is a single richness column)
    """

    if abundance_df.isna().sum().sum() > 0:
        richness_df = pd.DataFrame((~abundance_df.isna()).sum(axis=1))
    else:
        richness_df = pd.DataFrame((abundance_df != abundance_df.min().min()).sum(axis=1))
    richness_df.columns = ['richness']

    return richness_df


# microbiome strains
def add_strain_replacement(test, control=None, quantile=0.05):
    """
    Adds "replacement column that represents strain replacement to a dissimilarity data frame

    :param test: (pd.DataFrame) dissimilarity data frame with 'SampleName1/2' and 'dissimilarity' columns to which
    replacement column will be add to comparisons within individual
    :param control: (pd.DataFrame) dissimilarity data frame with 'SampleName1/2' and 'dissimilarity' columns from
    which strain replacement threshold will be calculated based on comparisons between individuals
    if not defined will use test data frame between individuals dissimilarity
    :param quantile: between individuals dissimilarity threshold from which replacement will be defined

    :return: test (pd.DataFrame)
    """

    # control data
    if control is None:
        control = test

    # otherwise changes inplace
    test = test.copy()
    control = control.copy()

    # remove duplicate comparisons and comparisons to the same sample
    control = control.reset_index()
    control = control[control['SampleName1'] != control['SampleName2']]
    control['s1'] = np.minimum(control['SampleName1'], control['SampleName2'])
    control['s2'] = np.maximum(control['SampleName1'], control['SampleName2'])
    control = control.drop_duplicates(['Species', 's1', 's2'])

    # remove within person comparisons
    if 'person1' in control.columns and 'person2' in control.columns:
        control = control[control['person1'] != control['person2']]
    else:
        print('the control samples are not identified by person, assuming all are between individuals comparisons')

    # get quantile
    control = control.groupby('Species')['dissimilarity'].quantile(quantile)

    # test data
    # remove comparisons to the same sample
    condition = (test.index.get_level_values('SampleName1') != test.index.get_level_values('SampleName2'))

    # remove between person comparisons
    if 'person1' in test.index.names and 'person2' in test.index.names:
        condition = condition & (test.index.get_level_values('person1') == test.index.get_level_values('person2'))
    else:
        print('the test samples are not identified by person, assuming all are within individual comparisons')

    test.loc[condition, 'replacement'] = test[condition].loc[control.index].groupby('Species')[
        'dissimilarity'].apply(lambda specie_data: specie_data > control.loc[specie_data.name])
    # control.index are control species

    return test


def add_unshared_pos(df):
    """
    Adds "unshared_pos1/2" columns of the unique avilable positions in each sample

    :param df: (pd.DataFrame) dissimilarity data frame with 'SampleName1/2' and 'shared_pos' columns

    :return: df (pd.DataFrame)
    """

    same_sample = df[df.index.get_level_values('SampleName1') == df.index.get_level_values('SampleName2')]
    same_sample = same_sample.reset_index().set_index(['Species', 'SampleName1'])[['shared_pos']].\
                              rename(columns={'SampleName1': 'SampleName', 'shared_pos': 'unshared_pos'})

    df = df.join(same_sample, on=['Species', 'SampleName1']).join(same_sample, on=['Species', 'SampleName2'],
            lsuffix='1', rsuffix='2')

    df['unshared_pos1'] = df['unshared_pos1'] - df['shared_pos']
    df['unshared_pos2'] = df['unshared_pos2'] - df['shared_pos']

    return df


# others

def segal_name(reps):
    """
    Get bacterias pretty name from Rep ID

    :param reps: (string or int or list) bacterias Rep ID, 'Rep_123' or 123, or list of those

    :return: (string) bacterias pretty name, 's. name'
    """

    if type(reps) == int or type(reps) == str:
        reps = [reps]

    for i in np.arange(len(reps)):
        reps[i] = int(str(reps[i]).replace('Rep_', ''))

    pretty_names = list(segal_df.loc[reps, 'pretty_name'].values)

    return pretty_names


def segata_name(SGB_ID):
    """
    Get bacterias full name from SGB ID
    
    :param SGB_ID: (string or int) bacterias SGB ID, 'SGB_123' or 123
    
    :return: (string) bacterias full name, 'k__?|p__?|c__?|o__?|f__?|g__?|s__?'
    """
    
    SGB_ID = int(str(SGB_ID).replace('SGB_', ''))

    full_name = segata_df.loc[segata_df['SGB ID'] == SGB_ID, 'Estimated taxonomy'].iloc[0].split('|')[-1]

    return full_name


def cat2binary(y):

    le = LabelEncoder()
    y = le.fit_transform(y)
    dict = {l: i for i, l in enumerate(le.classes_)}
    return y, dict


def mantel_test(s1, s2, s1_dis=True, s2_dis=False, maximal_filling=0.25, minimal_samples=20,
                method='pearson', permutations=100, alternative='two-sided'):
    # TODO: double check not sure logic make sense
    """
    Calculates the mantel test between two dissimilarity matrices

    :param s1: (pd.Series) first series values
    :param s2: (pd.Series) second series values
    :param s1_dis: (bool) whether the first series col values are already dissimilarity values
    :param s2_dis: (bool) whether the second series col values are already dissimilarity values
    :param maximal_filling: (float) maximal percentage of samples with missing values to fill with median value
    :param minimal_samples: (int) minimal number of samples to have in order to calculate mantel
    :param method: (str) correlation method 'pearson' or 'spearman'
    :param permutations: (int) number of permutations to do for the p_value
    :param alternative: (str) alternative hypothesis to use when calculating statistical significance:
    'two-sided', 'greater' or 'less'

    :return: r, p, n
    """

    def mantel_with_nans(x, y, method, permutations, alternative):

        """Compute correlation between distance matrices using the Mantel test.

        The Mantel test compares two distance matrices by computing the correlation
        between the distances in the lower (or upper) triangular portions of the
        symmetric distance matrices. Correlation can be computed using Pearson's
        product-moment correlation coefficient or Spearman's rank correlation
        coefficient."""

        # This is the skbio.stats.distance.mantel function with variations to accommodate missing values

        if method == 'pearson':
            corr_func = pearsonr
        elif method == 'spearman':
            corr_func = spearmanr
        elif method == 'kendalltau':
            corr_func = kendalltau
        else:
            raise ValueError("Invalid correlation method '%s'." % method)

        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")
        if alternative not in ('two-sided', 'greater', 'less'):
            raise ValueError("Invalid alternative hypothesis '%s'." % alternative)

        # x, y = _order_dms(x, y, strict=strict, lookup=lookup)

        n = x.shape[0]
        if n < 3:
            raise ValueError("Distance matrices must have at least 3 matching IDs "
                             "between them (i.e., minimum 3x3 in size).")

        x = DistanceMatrix(x.fillna(999))
        y = DistanceMatrix(y.fillna(999))
        x._data[x._data == 999] = np.nan
        y._data[y._data == 999] = np.nan

        x_flat = x.condensed_form()
        y_flat = y.condensed_form()

        na_loc = (np.isnan(x_flat)) | (np.isnan(y_flat))
        x_flat_no_na = x_flat[~na_loc]
        y_flat_no_na = y_flat[~na_loc]

        orig_stat = corr_func(x_flat_no_na, y_flat_no_na)[0]

        if permutations == 0 or np.isnan(orig_stat):
            p_value = np.nan
        else:
            permuted_stats = []
            for _ in range(permutations):
                # since the missing values shift relative places between the matrices,
                # different values will be dropped each time
                x_perm = x.permute(condensed=True)
                na_loc = (np.isnan(x_perm)) | (np.isnan(y_flat))
                x_flat_no_na = x_perm[~na_loc]
                y_flat_no_na = y_flat[~na_loc]
                permuted_stats.append(corr_func(x_flat_no_na, y_flat_no_na)[0])

            if alternative == 'two-sided':
                count_better = (np.absolute(permuted_stats) >=
                                np.absolute(orig_stat)).sum()
            elif alternative == 'greater':
                count_better = (permuted_stats >= orig_stat).sum()
            else:
                count_better = (permuted_stats <= orig_stat).sum()

            p_value = (count_better + 1) / (permutations + 1)

        return orig_stat, p_value, n

    def prep_series(s, s_dis):

        if s_dis:  # has dissimilarity values in a series form
            df = s.reset_index().pivot_table(index=[ind + '2' for ind in combined_indices],
                                             columns=[ind + '1' for ind in combined_indices], values='dissimilarity')
            df.index.names = combined_indices
            df.columns.names = combined_indices
            np.fill_diagonal(df.values, 0)
        else:  # calculate dissimilarity values from regular values
            ind2drop = list(set(s.index.names) - set(combined_indices))
            s.index = s.index.droplevel(ind2drop)
            df = pd.DataFrame(distance_matrix(x=np.reshape(s.values, (-1, 1)), y=np.reshape(s.values, (-1, 1))),
                              index=s.index, columns=s.index)

        df.index = df.index.reorder_levels(combined_indices)
        df.columns = df.columns.reorder_levels(combined_indices)

        # limit to samples that have at least maximal_filling percentage of existing dissimilarity measurements
        samples_mask = df.columns[df.isna().sum() < df.shape[0] * maximal_filling].values
        df = df.loc[samples_mask, samples_mask]

        # # find the dissimilarities that are still missing and fill them with medians
        # if df.isna().sum().sum() != 0:  # the condition is just for speed
        #     row_med = df.median(axis=0)
        #     col_med = df.median(axis=1)
        #     mean_medians = row_med.apply(lambda r_m: (r_m + col_med) / 2)
        #     df = df.mask(df.isna(), mean_medians)

        return df

    # sync indices names
    s1_ind = [ind[:-1] for ind in s1.index.names] if s1_dis else s1.index.names
    s2_ind = [ind[:-1] for ind in s2.index.names] if s2_dis else s2.index.names
    combined_indices = list(set(s1_ind) & set(s2_ind))

    # prep series to be a full dissimilarity series
    df1 = prep_series(s1, s1_dis)
    df2 = prep_series(s2, s2_dis)

    # sync indices values
    combined_indices = df1.index.intersection(df2.index)
    df1 = df1.loc[combined_indices, combined_indices]
    df2 = df2.loc[combined_indices, combined_indices]

    # mantel test
    if len(combined_indices) > minimal_samples:
        r, p, n = mantel_with_nans(df1, df2, method, permutations, alternative)
    else:
        r, p, n = None, None, None

    return r, p, n


def permutation_test(x1, x2, metric=np.average, n_permutations=100):
    """
    Calculate permutation test p_value for any metric

    :param x1: (any type of array) first data group
    :param x2: (any type of array) second data group
    :param metric: (func) any distribution metric like mean, median, kurtosis, etc. that returns a single value
    :param n_permutations: (int) number of permutations to perform

    :return: (float) p_val
    """

    # based on -
    # https://towardsdatascience.com/how-to-assess-statistical-significance-in-your-data-with-permutation-tests-8bb925b2113d
    # it also has correlation implementationon if needed

    x1 = np.array(x1)  # just to prevent problems
    x2 = np.array(x2)  # just to prevent problems

    # Bootstrapping - randomly sample without replacement two distributions with the size equal to the original
    # distributions from this pooled distribution to compute the absolute difference of the metric between the
    # two permuted samples. Repeated n_permutations times

    # Ground truth absolute difference between the labels from the two variables
    ground_truth = np.abs(metric(x1) - metric(x2))

    pooled_variables = list(x1) + list(x2)  # Pooled variables distribution
    pooled_shuffle = copy.copy(pooled_variables)  # Copy of the pooled variables distribution

    pooled_distribution = []  # Initialize permutation
    for i in np.arange(n_permutations):  # Permutation
        random.shuffle(pooled_shuffle)  # Shuffle the data
        # Permuted absolute difference of the two sampled distributions
        pooled_distribution.append(np.abs(metric(pooled_shuffle[:len(x1)]) - metric(pooled_shuffle[-len(x2):])))

    # the proportion of permuted differences higher than the ground truth difference is the significance value
    p_val = len(np.where(pooled_distribution >= ground_truth)[0]) / n_permutations

    return p_val


# Helper classes

class _Directories:

    # Initialization
    def __init__(self, base_directory):

        self.figs = os.path.join(base_directory, 'figs')
        self.excels = os.path.join(base_directory, 'excels')
        self.data_frames = os.path.join(base_directory, 'data_frames')
        self.models = os.path.join(base_directory, 'models')
        self.jobs = os.path.join(base_directory, 'jobs')

        directories = [self.figs, self.excels, self.data_frames, self.models, self.jobs]

        for Dir in directories:
            if not os.path.exists(Dir):
                os.mkdir(Dir)


class _Parameters:

    # Initialization
    def __init__(self, study=None, controls=None, colors=None,
                 alpha=None, detection_threshold=None, dissimilarity_threshold=None):

        self.study = study
        self.controls = controls if controls is not None else {}
        self.colors = colors if colors is not None else {}
        self.alpha = alpha
        self.detection_threshold = detection_threshold
        self.dissimilarity_threshold = dissimilarity_threshold

        # in case some key parameters were not set define them with default values
        defaults = ['group', 'time_point']
        for param in defaults:
            if param not in self.controls.keys():
                self.controls[param] = 'fake_{}'.format(param)
            if param not in self.colors.keys():
                self.colors['fake_{}'.format(param)] = 'gray'
        self.colors[''] = 'gray'  # necessary


if __name__ == '__main__':
    pass
