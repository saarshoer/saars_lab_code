# Imports
# files
import os
import pickle
import tarfile
from glob import glob
from ftplib import FTP

# data
import math
import numpy as np
import pandas as pd

# statistics
from mne.stats.multi_comp import fdr_correction
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_rel, ttest_1samp, binom_test, spearmanr

# models
from xgboost import XGBClassifier, XGBRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold, train_test_split, GridSearchCV

# dimensionality reduction
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# clustering
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster

# plots
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm

# lab
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders import GutMBLoader, OralMBLoader, BloodTestsLoader, BodyMeasuresLoader, \
    CGMLoader, DietLoggingLoader

segata_df = pd.read_csv(
            '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Segata/SupplementaryTable8-SGBsDescription.csv')

# Main class
class Study:

    # Initialization
    def __init__(self,

                 # Default parameters
                 study=None, controls=None, colors=None,
                 alpha=0.05, detection_threshold=0.0001, dissimilarity_threshold=1/20000,
                 base_directory=os.getcwd(), figsize_factor=1.5):

        # Parameters
        self.params = _Parameters(study=study, controls=controls, colors=colors, figsize_factor=figsize_factor,
                                  alpha=alpha,
                                  detection_threshold=detection_threshold,
                                  dissimilarity_threshold=dissimilarity_threshold)

        # Directories
        self.dirs = _Directories(base_directory)

        # Objects
        self.objs = {}

        # Figures
        rcParams['figure.figsize'] = (rcParams['figure.figsize'][0] * figsize_factor,
                                      rcParams['figure.figsize'][1] * figsize_factor)

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
        :param columns_from_file: (str, list or dict) to add from
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
        if type(data) == str:
            if data == 'gutMB':
                lab_data = GutMBLoader.GutMBLoader().get_data('segata_species', study_ids=self.params.study)
                fig_abundance_reads(lab_data.df_metadata)
            elif data == 'oralMB':
                lab_data = OralMBLoader.OralMBLoader().get_data('segata_species', study_ids=self.params.study)
                fig_abundance_reads(lab_data.df_metadata)
            elif data == 'blood':
                lab_data = BloodTestsLoader.BloodTestsLoader().get_data(study_ids=self.params.study)
            elif data == 'body':
                lab_data = BodyMeasuresLoader.BodyMeasuresLoader().get_data(study_ids=self.params.study)
            elif data == 'cgm':
                lab_data = CGMLoader.CGMLoader().get_data(study_ids=self.params.study)
            elif data == 'diet':
                loader = DietLoggingLoader.DietLoggingLoader()
                lab_data = loader.get_data(study_ids=self.params.study)
                lab_data.df = loader.add_nutrients(lab_data.df,
                                                   nutrient_list=['protein_g', 'totallipid_g', 'carbohydrate_g'])
            else:
                raise Exception('data string is not valid')

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
                # read file to df according to the file's type
                file_extension = os.path.splitext(file)[1]
                if file_extension == '.xlsx':
                    file_df = pd.read_excel(file)
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
                file_df['RegistrationCode'] = file_df['RegistrationCode'].astype(str)

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

        # index manipulations
        df = df.set_index(added_columns, append=True)
        df = df.rename_axis(index=default_name_changes)

        indices_order = list(df.index.names)
        indices_order.sort(key=lambda i: (default_indices_order + list(df.index.names)).index(i))
        df = df.reorder_levels(indices_order)
        df = df.sort_values(by=df.index.names)

        return df

    def prep_df(self, obj, time_points=None, mean_time_point=False, n_entries=None,
                indices_dict=None, log=True):
        """
        Handel the time_point, replace indices values, handel missing values (if abundance) and optionally log10

        :param obj:
        :param time_points: (list) to include
        :param mean_time_point: (bool) if to mean between duplicate entries of the same time_point
        :param n_entries: (int) amount of entries to accept per person
        :param indices_dict: (dict of dicts) outer key for index name
        inner key for index value and value after replacement
        :param log: (bool) whether or not to log10 the values

        :return: (pd.DataFrame) df
        """

        def fig_abundance_distribution(abundance_df):

            df1 = pd.DataFrame(abundance_df.stack()).reset_index()
            df2 = pd.DataFrame(np.log10(abundance_df).stack()).reset_index()
            df1['log'] = False
            df2['log'] = True

            # if there is no actual control time_point we have a problem
            if self.params.controls['time_point'] != 'fake_time_point':
                df3 = pd.DataFrame(self.get_delta_df(abundance_df, self.params.controls['time_point'])
                                   .stack()).reset_index()
                df4 = pd.DataFrame(self.get_delta_df(np.log10(abundance_df), self.params.controls['time_point'])
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
                                                .stack()).reset_index().rename(columns={0: 'n_species'})

            if self.params.controls['time_point'] == 'fake_time_point':
                n_species_per_sample['time_point'] = 'fake_time_point'

            # hue 'group'
            g = sns.FacetGrid(n_species_per_sample, col='time_point', hue='group', palette=self.params.colors,
                              sharex=True, sharey=True, margin_titles=True)
            g = g.map(sns.distplot, 'n_species', hist=True, kde=False, bins=10)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            g.set_ylabels('n_samples')
            title = '{}species distribution colored by group'.format(obj.type.replace('abundance', ''))
            plt.suptitle(title)
            plt.subplots_adjust(top=0.8)

            plt.savefig(os.path.join(self.dirs.figs, title))

            # hue 'time_point'
            g = sns.FacetGrid(n_species_per_sample, col='group', hue='time_point', palette=self.params.colors,
                              sharex=True, sharey=True, margin_titles=True)
            g = g.map(sns.distplot, 'n_species', hist=True, kde=False, bins=10)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            g.set_ylabels('n_samples')
            title = '{}species distribution colored by time_point'.format(obj.type.replace('abundance', ''))
            plt.suptitle(title)
            plt.subplots_adjust(top=0.8)

            plt.savefig(os.path.join(self.dirs.figs, title))

        def time2time_point(person_abundance_df):

            # add time_point by the time of the points, for example: 03/05, 01/04 will become 1, 0
            person_abundance_df['time_point'] = person_abundance_df.index.get_level_values('time').argsort()

            return person_abundance_df

        def declare_missing_values(person_abundance_df):

            # in places where all column values equal to detection threshold declare them as missing values
            person_abundance_df.loc[:, (person_abundance_df == self.params.detection_threshold).all()] = np.nan

            return person_abundance_df

        # retrieving
        df = obj.df

        # add time_point by the time of the points, for example: 03/05, 01/04 will become 1, 0
        if 'time_point' not in df.index.names:
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
            # this obliviates time index because the mean function can not handel datetime values

        # filter out cases where a person has the wrong number of samples
        if n_entries is not None:
            df = df.groupby('person').filter(lambda x: len(x) == n_entries)

        # replace indices values according to dictionary
        if indices_dict is not None:
            for ind in indices_dict.keys():
                if ind in df.index.names:
                    df = df.rename(index=indices_dict[ind], level=ind)

        if 'abundance' in obj.type:
            # fill missing values caused by MBLoader with detection threshold
            df = df.fillna(self.params.detection_threshold)

            # filter out empty species before declare_missing_values to shorten run time
            df = df.loc[:, (df != self.params.detection_threshold).any()]

            # declare missing values
            if len(np.unique(df.index.get_level_values('time_point'))) > 1:
                df = df.groupby(['person']).apply(declare_missing_values)
            else:
                df = df.replace(self.params.detection_threshold, np.nan)

            # figures
            fig_abundance_distribution(df)
            fig_species_distribution(df)

        # filter out empty columns (species/test)
        df = df.dropna(axis=1, how='all')

        # convert the values to the order of magnitude
        if log:
            df = np.log10(df)

        return df

    # analysis
    def comp_stats(self, obj, test, between, delta=False, minimal_samples=0.1, internal_use=False):
        """
        Compute the statistical significance for the difference between elements

        :param obj: (object) object to calculate the statistics on
        :param test: (str) mannwhitneyu, wilcoxon, ttest_ind, ttest_rel, ttest_1samp, binom_test
        :param between: (str) group or time_point
        :param delta: (bool) whether or not to calculate the statistics based on the difference between
        time_point values
        :param minimal_samples: (percentage or int)
        :param internal_use: (bool) whether an internal function is running it

        :return: None
        """

        # TODO: maybe in case of abundance remove detection threshold form comparisons that do not include delta or paired

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
                                     sharey=True)  # sharex

                plt.subplots_adjust(bottom=0.5, left=0.14)
                plt.suptitle('{}\n{} - significant results'.format(obj.type, test))

            ax = axes_internal[major_elements.index(major_e) * len(minor_elements) + minor_elements.index(minor_e)]

            # in case there are some significant results
            curr_data_df = curr_data_df.dropna()
            if curr_data_df.shape[0] != 0:

                # transfer the columns to rows with identifiers
                # reset index behaves differently for single level versus multilevel indices
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
                ax.set_ylabel(y_label)
                ax.set_xlabel('')
                ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=90)
                ax.legend().set_visible(False)
                if delta:
                    xmin, xmax = ax.get_xlim()
                    ax.hlines(y=0, xmin=xmin, xmax=xmax, colors='lightgray')

            else:
                # texting
                ax.text(x=(ax.get_xlim()[0]+ax.get_xlim()[1])/2, y=(ax.get_ylim()[0]+ax.get_ylim()[1])/2,
                        s='no significant results', horizontalalignment='center')
                # TODO: fix bug that if this is the first ax then the middle of the graph is changed when plotting the second ax
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
                figure_internal = None
                axes_internal = None

            return figure_internal, axes_internal

        # retrieving the data frames from the objects
        df = obj.df

        # change the values in the data frames to be the change in value between two time points:
        # any time point and the control time point
        if delta:
            df = self.get_delta_df(df, self.params.controls['time_point'])
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
        if major == 'group':
            major_control = self.params.controls['group']
        elif major == 'time_point':
            major_control = self.params.controls['time_point']
        else:  # just to prevent a warning later on
            major_control = None

        # tests without control data
        tests_without_control = ['ttest_1samp', 'binom_test']

        # data frame to hold all the data compered for further (figure) analysis
        columns = pd.MultiIndex.from_product([major_elements, minor_elements], names=[major, minor])
        data_df = pd.DataFrame(index=df.columns, columns=columns)
        figure = None
        axes = None

        # remove control from the list of elements
        if test not in tests_without_control:
            if len(major_elements) < 2:
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
                        # mediterranean handling just to limit excel sheet name to the 31 character count limit
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

                    # define the test data not matter the test
                    test_data = df[col] \
                        .xs(major_e, level=major, drop_level=False) \
                        .xs(minor_e, level=minor, drop_level=False).dropna()
                    data_df.loc[col, (major_e, minor_e)] = test_data.values

                    # handling cases where there is missing value in one time_point but not in the other
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

                # fdr correction
                _, stats_df_list[-1]['p_FDR'] = fdr_correction(stats_df_list[-1]['p'], alpha=self.params.alpha)
                stats_df_list[-1].sort_values('p_FDR', inplace=True)

                # delete from the data_df un-significant results
                significant_indices = stats_df_list[-1][stats_df_list[-1]['p_FDR'] < self.params.alpha].index
                non_significant_indices = data_df.index.difference(significant_indices)
                data_df.loc[non_significant_indices, (major_e, minor_e)] = np.nan

                # print for the user
                print('')
                print(test)
                print(stats_df_list[-1].name)
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
                                                  curr_data_df=data_df.xs(major_e, level=major, axis=1, drop_level=False)
                                                  .xs(minor_e, level=minor, axis=1))

        if not internal_use:
            # excel finishes
            excel_writer.save()
            excel_writer.close()
        else:
            return stats_df_list

    def score_models(self, xobj, yobj, delta=True, minimal_samples=0.9,
                     model_type='linear', n_repeats=1, n_splits=5, random_state=None, hyper_parameters=None,
                     send2queue=False, save=False):
        """
        Calculate n_iterations prediction models or grid search for the hyper parameters for the yobj based on the xobj
        and return the summary statistics on the models scores

        :param xobj: (object) the features object
        :param yobj: (object) the target object
        :param delta: (bool) whether or not to calculate the models based on the difference between time_point values
        :param minimal_samples: (percentage or int) that needs to be included in the models
        and therefore cannot be deleted in the missing values filtering process
        :param model_type: (str) 'linear' or 'xgb'
        :param n_repeats: (int) number of models to be created for each target
        :param n_splits: (int) number of time to split the data in the KFold
        :param random_state: (int) seed for the random process used for splitting the data
        :param hyper_parameters: (dict) dictionary with hyper parameters for the model {key: [value/s]}
        if value is single will be set directly, if values are multiple will perform a grid search
        :param send2queue: (bool) whether or not to send the jobs to the queue or to run locally
        :param save: (bool) whether or not to save the models

        :return: (pd.DataFrame) scores_df
        """

        # TODO: add SHAP analysis
        def score_model(x, y):

            # check if column is classifiable in order to know which model to use
            if sorted(np.unique(y)) == [0, 1] or sorted(np.unique(y)) == [False, True] \
                    or not all(np.issubdtype(val, np.number) for val in y):
                classifiable = True
            else:
                classifiable = False

            # skip column in some cases
            if model_type == 'linear' and classifiable:
                # print('skipping classifiable column while model is linear - {}'.format(col))
                return [np.nan, np.nan]
            if any(np.unique(y) == np.nan):
                # print('skipping column with {} missing values - {}'.format(y_df[col].isna().sum(), col))
                return [np.nan, np.nan]
            if np.unique(y).shape[0] == 1:
                # print('skipping column with a single value ({}) - {}'.format(y_df[col].unique(), col))
                return [np.nan, np.nan]

            # a list to contain the iterations scores
            scores = []

            # empty array to hold all predictions in order to score them later together
            # this will be overwritten in case of a grid search
            y_pred = np.empty(y.shape)
            y_pred.fill(np.nan)

            # linear model
            if model_type == 'linear':
                model_instance = LinearRegression(**direct_hyper_params)
                score = 'r2'

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
            kf = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)  # shuffle=True

            # grid search mode
            if len(search_hyper_params) != 0:

                model = GridSearchCV(model_instance, search_hyper_params, scoring=score, cv=kf)

                x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.5, random_state=random_state)

                model.fit(x_train, y_train)

                scores.append(model.score(x_test, y_test))

            # not grid search mode
            else:

                model = model_instance

                for train_index, test_index in kf.split(x):
                    x_train, x_test = x[train_index], x[test_index]
                    y_train, y_test = y[train_index], y[test_index]

                    model.fit(x_train, y_train)

                    scores.append(model.score(x_test, y_test))

            # save/load the model
            if save:
                pickle.dump(model, open(os.path.join(self.dirs.models,
                            'model {} {} {}{} {}.sav'.format(xobj.type, yobj.type, col, delta_str, model_type)), 'wb'))
                # model = pickle.load(open('model.sav', 'rb'))

            return [np.mean(scores), np.std(scores)]

        # retrieving the data frames from the objects
        if type(xobj) == self.Object and type(yobj) == self.Object:
            x_df = xobj.df
            y_df = yobj.df
        else:
            raise Exception('xobj and yobj type is not object')

        # change the values in the data frames to be the change in value between two time points:
        # any time point and the control time point
        if delta:
            x_df = self.get_delta_df(x_df, self.params.controls['time_point'])
            y_df = self.get_delta_df(y_df, self.params.controls['time_point'])
            delta_str = ' delta'  # for file names
        else:
            delta_str = ''  # for file names

        # convert the minimal_samples percentage from the argument to the number of samples
        if minimal_samples < 1:  # meaning minimal_samples is in percentage
            minimal_samples = round(minimal_samples * x_df.shape[0])
            # needs to be after the delta because delta effects the shape

        # deal with missing values
        if 'abundance' in xobj.type:  # only for abundance otherwise the min min might not be the best solution
            x_df = x_df.fillna(x_df.min().min())
        # drop samples with missing values as long as it does not go over the minimal_samples
        na_columns = (y_df.isna().sum(axis=0) > 0) & (y_df.isna().sum(axis=0) < (x_df.shape[0] - minimal_samples))
        # columns with small number of missing values and thus the entries missing in them can be deleted
        na_person = y_df.loc[:, na_columns].isna().groupby('person').any().any(axis=1)
        person2delete = na_person[na_person.values].index.values
        # person from rows that can be deleted by the previous logic
        y_df = y_df[~y_df.index.isin(person2delete, level='person')]

        # find the indices that exists in both the x_df and the y_df
        combined_indices = x_df.index.get_level_values('person').intersection(y_df.index.get_level_values('person'))
        # filter the data frames to only include these indices
        x_df = x_df[x_df.index.isin(combined_indices, level='person')]
        y_df = y_df[y_df.index.isin(combined_indices, level='person')]

        # conversion
        x_ = np.array(x_df)

        # split hyper parameters to those given and those that need to be optimized
        if hyper_parameters is None:
            hyper_parameters = {}
        # split the hyper parameters to ones you can directly give the model instance and to ones you need to search for
        direct_hyper_params = {key: value[0] for key, value in hyper_parameters.items() if len(value) == 1}
        search_hyper_params = {key: value for key, value in hyper_parameters.items() if len(value) != 1}

        # random_state is used by the xgb model to split the data the same way in the KFold
        # if random_state is not None:
        #     np.random.seed(random_state)
        # random_state = np.random.randint(100, size=n_iterations)  # the 100 is arbitrary

        # queue
        original_dir = os.getcwd()
        os.chdir(os.path.join(self.dirs.excels, '..', 'jobs'))
        queue = qp if send2queue else fakeqp

        with queue(jobname='model', _delete_csh_withnoerr=True, q=['himem7.q']) as q:

            q.startpermanentrun()
            tkttores = {}

            # create a model for each target column
            for col in y_df.columns:

                if y_df[col].isna().sum() != 0:
                    print('skipping column with {} missing values - {}'.format(y_df[col].isna().sum(), col))
                    continue
                if y_df[col].unique().shape[0] == 1:
                    print('skipping column with a single value ({}) - {}'.format(y_df[col].unique(), col))
                    continue

                # conversion
                y_ = np.array(y_df[col])

                # send job
                tkttores[col] = q.method(score_model, (x_, y_))

            # summarize the scores using the iterations scores distribution if exist
            scores_df = pd.DataFrame(index=y_df.columns, columns=['mean', 'std'])

            for k, v in tkttores.items():
                scores_df.loc[k, ['mean', 'std']] = q.waitforresult(v)

        os.chdir(original_dir)

        # drop all the un-ran columns
        scores_df = scores_df.dropna(axis=0, how='all')

        # sort values by score
        scores_df = scores_df.sort_values('mean', ascending=False)

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
        The resulting excel is filtering while the figures are not

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

                    # len 2 to be plottable, shape > 0 to have values, param different than hue for hue to make sense
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
            x_df = self.get_delta_df(x_df, self.params.controls['time_point'])
            y_df = self.get_delta_df(y_df, self.params.controls['time_point'])
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

    def dim_reduction(self, obj, n_pca_comp=50, n_tsne_comp=2):
        """
        Reduces the dimensionality of a features X samples data frame
        and save the best 2 components as figures, colored by each index

        :param obj: (Object) with data frame containing features X samples,
        index levels are used to color and create different figures
        :param n_pca_comp: (None, 0 or int) number of desired components out of the pca
        :param n_tsne_comp: (None, 0 or int) number of desired components out of the pca

        ** not filling either of the n_*_comp will use the dim_reduction default arguments
        ** filling either of the n_comp with 'None' will use pca/tsne default arguments
        ** filling either of the n_comp with '0' will cause the function to skip the method (pca/tsne)

        :return: None
        """

        # the default n_pca_comp is 50 because that is the maximal number of features recommend to the tsne function
        # the default n_tsne_comp is 2 because this is plotable

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
                if n_pca_comp != 0:
                    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1],
                                    hue=hue, palette=palette, alpha=0.3, ax=axes[0])
                    axes[0].set_title('pca')
                    axes[0].set_xlabel('pca component 1 - {:.2%}'.format(pca.explained_variance_ratio_[0]))
                    axes[0].set_ylabel('pca component 2 - {:.2%}'.format(pca.explained_variance_ratio_[1]))

                # tsne
                if n_tsne_comp != 0:
                    sns.scatterplot(x=tsne_result[:, 0], y=tsne_result[:, 1],
                                    hue=hue, palette=palette, alpha=0.3, ax=axes[1])
                    axes[1].set_title('tsne')
                    axes[1].set_xlabel('tsne component 1')
                    axes[1].set_ylabel('tsne component 2')

                title = 'dim reduction {} colored by {}'.format(obj.type, index_name)
                fig.suptitle(title)
                plt.savefig(os.path.join(self.dirs.figs, title))

        # retrieving the data frames from the object
        if type(obj) == self.Object:
            df = obj.df
        else:
            raise Exception('obj is not object')

        if n_pca_comp == 0 and n_tsne_comp == 0:
            raise Exception('both reduction methods are empty')

        # pca
        if n_pca_comp != 0:
            pca = PCA(n_components=n_pca_comp)
            pca_result = pca.fit_transform(df.values)
        else:
            pca_result = df.values

        # tsne
        if n_tsne_comp != 0:
            tsne = TSNE(n_components=n_tsne_comp)
            tsne_result = tsne.fit_transform(pca_result)
        else:
            tsne_result = pca_result

        fig_best_components()

    def fig_snp_heatmap(self, obj, maximal_filling=0.25, minimal_samples=20,
                        annotations=None, cmap=None, log_colors=False):
        """
        Plot a heatmap based on the SNP dissimilarity data frame for each species

        :param obj: (Object) with SNP dissimilarity data frame
        :param maximal_filling: (float) maximal percentage of samples with missing values to fill with median value
        JUST for clustering, it will not appear in the heatmap!!
        :param minimal_samples: (int) minimal number of samples to have in order to draw
        :param annotations: (list of str) columns to annotate by
        :param cmap: (str) dissimilarity color map name
        :param log_colors: (bool) whether to set the heatmap colors to log scale

        :return: None
        """

        if annotations is None:
            annotations = ['group', 'time_point']

        # create all the annotation data for the colored bars
        annotations_df = obj.df
        annotations_df = annotations_df.reset_index()[[s for s in annotations_df.index.names if '1' in s]] \
            .set_index('SampleName1').drop_duplicates()  # take only sample1 metadata
        annotations_df.columns = annotations_df.columns.str.rstrip('1')  # remove 1 suffix from the names
        annotations_df = annotations_df[annotations]  # limit the annotation to these categories
        annotations_df = annotations_df.iloc[~annotations_df.index.duplicated()]  # drop duplicated samples
        annotations_labels = np.unique(annotations_df.values)
        colors_df = annotations_df.replace(self.params.colors)  # replace values with corresponding colors

        # for each species
        for species in obj.df.index.get_level_values('Species').unique():

            # create a square data matrix (drop metadata)
            df = obj.df.loc[species].reset_index().pivot_table(index='SampleName2',
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
                                         method='average', optimal_ordering=True)
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
                    g = sns.clustermap(df, mask=(df >= 0),  # the mask here is for everything
                                       xticklabels=False, yticklabels=False,
                                       row_linkage=df_linkage, col_linkage=df_linkage,
                                       row_colors=colors_df, col_colors=colors_df)

                    norm = LogNorm(g.data2d.min().min(), g.data2d.max().max(), clip=False) if log_colors else None

                    sns.heatmap(g.data2d, mask=na_mask.loc[g.data2d.index, g.data2d.columns],
                                xticklabels=False, yticklabels=False,
                                ax=g.ax_heatmap, cbar_ax=g.cax,
                                cmap=cmap, norm=norm,
                                cbar_kws={'label': 'dissimilarity', 'orientation': 'horizontal', 'ticks': None})
                    # cbar ticks are None because otherwise in log scale there are no ticks

                    # statistical test
                    stats_df = annotations_df.loc[df.iloc[g.dendrogram_col.reordered_ind].index]
                    stats_df['rank'] = np.arange(stats_df.shape[0])

                    text = ''

                    for anno in annotations:

                        # statistics
                        if len(np.unique(stats_df[anno])) == 2:
                            s, p = mannwhitneyu(
                                x=stats_df.loc[stats_df[anno] == np.unique(stats_df[anno])[0], 'rank'].values,
                                y=stats_df.loc[stats_df[anno] == np.unique(stats_df[anno])[1], 'rank'].values,
                                use_continuity=True, alternative='two-sided')  # TODO: think about use_continuity
                            p = round(p, 3)
                        else:
                            p = 'NA'

                        # rand index between the annotation and the clusters
                        org_clusters = annotations_df.loc[df.index, anno].values.flatten()
                        new_clusters = fcluster(df_linkage, t=len(np.unique(org_clusters)),
                                                criterion='maxclust')

                        RI = round(adjusted_rand_score(org_clusters, new_clusters), 3)

                        # TODO: think what should be the t for inconsistent fcluster
                        # new_clusters2 = fcluster(df_linkage, t=1, criterion='inconsistent')
                        # rand_index2 = adjusted_rand_score(org_clusters, new_clusters2)

                        text = '{}\n{}\np={}, RI={}'.format(text, anno, p, RI)

                    g.fig.text(0.775, 0.875, text)
                    # TODO: find a better way to position the text

                    title = '{}\n{}'.format(obj.type, species)
                    g.fig.suptitle(title)

                    # annotations legend
                    for label in annotations_labels:
                        g.ax_col_dendrogram.bar(0, 0, color=self.params.colors[label], label=label, linewidth=0)

                    # position of annotations legend
                    g.ax_col_dendrogram.legend(ncol=1, frameon=False, bbox_to_anchor=(-0.035, 1))

                    # position of dissimilarity legend
                    g.cax.set_position([.35, .05, .5, .01])

                    if not os.path.exists(os.path.join(self.dirs.figs, obj.type)):
                        os.makedirs(os.path.join(self.dirs.figs, obj.type))
                    plt.savefig(os.path.join(self.dirs.figs, obj.type, title.replace('\n', ' ')), pad_inches=0.5)
                    plt.close()

    def fig_snp_scatter_box(self, obj, subplot='group', minimal_comparisons=45, species=25, height=12, aspect=0.5):
        """
        Plot a the dissimilarity distribution between and within people based on the SNP dissimilarity data frame
        for each species

        :param obj: (Object) with SNP dissimilarity data frame
        :param subplot: (str) group or time_point as subplots
        :param minimal_comparisons: (int) minimal number of comparisons between people per specie within a subplot
        :param species: (int or list) int of number of species per plot or a list of specific species to plot
        :param height: (int) sns.FacetGrid height argument
        :param aspect: (float) sns.FacetGrid aspect argument

        :return: None
        """

        def scatter_box_plot(**kwargs):
            # retract the data
            sub_data = kwargs.pop('data')

            # within individuals
            sns.scatterplot('dissimilarity', 'Species', major,
                            data=sub_data.loc[same_person],
                            hue_order=np.unique(sub_data.loc[same_person, major]), alpha=0.2, legend='brief', **kwargs)

            # between individuals
            sns.boxplot('dissimilarity', 'Species', minor,
                        data=sub_data.loc[~same_person & same_minor],
                        hue_order=np.unique(sub_data.loc[~same_person & same_minor, minor]),
                        width=0.5, fliersize=0, **kwargs)

            # between individuals - significance
            sns.scatterplot(1/15, 'Species', minor,
                            data=stats_df.loc[stats_df[major].isin(np.unique(sub_data[major])) &
                                              stats_df['Species'].isin(np.unique(sub_data['Species']))],
                            hue_order=np.unique(stats_df.loc[stats_df[major].isin(np.unique(sub_data[major])) &
                                                stats_df['Species'].isin(np.unique(sub_data['Species'])), minor]),
                            marker='*', s=120, legend=False, **kwargs)

        # data manipulation
        df = obj.df
        df = df.reset_index()
        df = df[df['SampleName1'] != df['SampleName2']]
        df['s1'] = np.minimum(df['SampleName1'], df['SampleName2'])
        df['s2'] = np.maximum(df['SampleName1'], df['SampleName2'])
        df = df.drop_duplicates(['Species', 's1', 's2'])

        # can subplot groups (control/test1/test2/...) or time points (0/1/2...)
        subplot_options = ['group', 'time_point']  # this list has to be length 2
        if subplot not in subplot_options:
            raise Exception('between not valid')

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

        # species argument
        if type(species) == list:  # by specific species
            df = df[df['Species'].isin(species)]
            figures = [0]
            species = len(species)

        # statistics
        obj4stats = Study.Object(obj_type='{}_'.format(obj.type), columns='bacteria')
        obj4stats.df = df.loc[~same_person & same_minor].set_index(['Species', major, minor], append=True)
        obj4stats.df = obj4stats.df['dissimilarity'].unstack('Species')

        stats_df = self.comp_stats(obj4stats, test='mannwhitneyu', between=minor,
                                   minimal_samples=minimal_comparisons, internal_use=True)

        for i in np.arange(len(stats_df)):
            stats_df[i][major] = stats_df[i].name.split('__')[1]
            stats_df[i][minor] = stats_df[i].name.split('__')[0]
        stats_df = pd.concat(stats_df)

        stats_df = stats_df[stats_df['p_FDR'] < self.params.alpha].reset_index()

        # removing species that lack inter or intra person data
        df = df[df['Species'].isin(list(
            set(df.loc[same_person, 'Species']).intersection(set(df.loc[~same_person & same_minor, 'Species']))))]

        df = df[df['Species'].isin(list(
            set(df['Species']).intersection(set(stats_df['Species']))))]
        stats_df = stats_df[stats_df['Species'].isin(list(
            set(df['Species']).intersection(set(stats_df['Species']))))]

        # species argument
        if type(species) == int:
            figures = np.arange(0, len(np.unique(df['Species'])), species)

        # split to multiple figures
        for i_figure in figures:

            sub_df = df[df['Species'].isin(list(df.groupby('Species').groups)[i_figure:i_figure+species])]

            # plotting
            g = sns.FacetGrid(sub_df, col=major, col_order=np.unique(df[major]), height=height, aspect=aspect)
            g = g.map_dataframe(scatter_box_plot, palette=self.params.colors)

            g.add_legend()
            plt.xscale('log')
            plt.xlim([self.params.detection_threshold/3, 10**-1])

            g.set_axis_labels(x_var='dissimilarity', y_var='species')
            g.set_titles(row_template='{row_name}', col_template='{col_name}')
            for ax in g.axes.flatten():
                if ax.get_title() in self.params.colors.keys():
                    color = self.params.colors[ax.get_title()]
                else:
                    color = 'black'
                ax.set_title(label=ax.get_title(), color=color)
            title = '{} distribution by {}'.format(obj.type, subplot)
            g.add_legend()
            plt.suptitle(title)
            plt.tight_layout()
            plt.subplots_adjust(top=0.88)

            plt.savefig(os.path.join(self.dirs.figs, '{} s{}-{}'.format(title, i_figure, i_figure+species)))

    # data frame computations

    @staticmethod
    def get_diversity_df(abundance_df):
        """
        Compute the Shanon's alpha diversity index (using log10) based on an "abundance_df"

        :param abundance_df: (pd.DataFrame) data frame to compute on, if not given takes the self.objs.abundance.df

        :return: (pd.DataFrame) "diversity_df"
        (same as abundance_df just instead of bacteria columns there is a single diversity column)
        """

        # Shannon's alpha diversity index = -sum(Pi*log10(Pi))

        # revert values to their pre log state
        if not (0 <= abundance_df.min().min() and abundance_df.max().max() <= 1):
            abundance_df = (10 ** abundance_df)

        diversity_df = pd.DataFrame(-(abundance_df * np.log2(abundance_df)).sum(axis=1).dropna())
        diversity_df.columns = ['diversity']

        return diversity_df

    @staticmethod
    def get_delta_df(regular_df, control_time_point):
        """
        Subtract from each time point values the self.params.control_time values

        :param regular_df: (pd.DataFrame) data frame to calculate for the delta
        :param control_time_point: time point values to subtract for all other time points values

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
            delta_df.loc[idx] = \
                (regular_df.xs(time_point, level='time_point') -
                 regular_df.xs(control_time_point, level='time_point')).values
            delta_df = delta_df.rename(
                index={time_point: '{}-{}'.format(time_point, control_time_point)},
                level='time_point')

        # delete the control time points
        idx = delta_df.xs(control_time_point, level='time_point', drop_level=False).index
        delta_df = delta_df.drop(idx)

        return delta_df

    @staticmethod
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

    @staticmethod
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


def bac_full_name(SGB_ID):
    """
    Get bacterias full name from SGB ID
    
    :param SGB_ID: (string or int) bacterias SGB ID, 'SGB_123' or 123
    
    :return: (string) bacterias full name, 'k__?|p__?|c__?|o__?|f__?|g__?|s__?'
    """
    
    SGB_ID = int(str(SGB_ID).replace('SGB_', ''))

    full_name = segata_df.loc[segata_df['SGB ID'] == SGB_ID, 'Estimated taxonomy'].iloc[0]

    return 'SGB_{} {}'.format(SGB_ID, full_name)


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
    def __init__(self, study=None, controls=None, colors=None, figsize_factor=None,
                 alpha=None, detection_threshold=None, dissimilarity_threshold=None):

        self.study = study
        self.controls = controls
        self.colors = colors
        self.figsize_factor = figsize_factor
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


if __name__ == "__main__":

    print(help(Study))
