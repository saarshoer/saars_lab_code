# Imports
# files
import os
import tarfile
from glob import glob
from ftplib import FTP

# data
import numpy as np
import pandas as pd

# models
from sklearn.model_selection import KFold
from xgboost import XGBClassifier, XGBRegressor
from sklearn.linear_model import LinearRegression

# statistics
from mne.stats.multi_comp import fdr_correction
from sklearn.metrics import roc_auc_score, r2_score
from scipy.stats import mannwhitneyu, ttest_ind, ttest_rel, ttest_1samp, binom_test, pearsonr, spearmanr

# plots
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# lab
from LabData.DataLoaders import MBLoader, BloodTestsLoader, BodyMeasuresLoader, CGMLoader

rcParams['figure.figsize'] = (rcParams['figure.figsize'][0]*1.5, rcParams['figure.figsize'][1]*1.5)


# Main class
class Study:

    # Initialization
    def __init__(self,

                 # Default parameters
                 study=None, controls=None, colors=None,
                 alpha=0.05, detection_threshold=0.0001,
                 base_directory=os.getcwd()):

        # Parameters
        self.params = _Parameters(study=study, controls=controls, colors=colors,
                                  alpha=alpha, detection_threshold=detection_threshold)

        # Directories
        self.dirs = _Directories(base_directory)

        # Objects
        self.objs = _Objects()

    # Functions

    # data generation
    def load_df(self, data,
                columns_from_metadata=None,
                columns_from_file=None, file=None, file_index=None):
        """
        Load a data frame and add columns from external sources

        :param data: (str or pd.DataFrame) string of data type ('microbiome', 'blood', 'body' or 'cgm') or a data frame
        :param columns_from_metadata: (str, list or dict) to add from the LabDataLoader meta data data frame
        :param columns_from_file: (str, list or dict) to add from
        :param file: (str or pd.DataFrame) from which to add columns ('xlsx' or 'csv')
        :param file_index: (str, list or dict) index in file to join by

        :return: (pd.DataFrame) df
        """

        def fig_abundance_reads(input_metadata_df):
            columns = ['RawRC', 'PostTrimRC', 'PostQCRC', 'PostHGFRC', 'UnalignedRC']

            plt.figure()
            sns.boxplot(data=input_metadata_df[columns])
            plt.title('Read Count')
            plt.ylabel('read count')
            plt.ylim([0, 0.4 * (10 ** 8)])

            plt.savefig(os.path.join(self.dirs.figs, 'read count'))

        # default parameters to always act on
        default_indices_order = ['group', 'person', 'time_point', 'time', 'sample']
        default_name_changes = {'SampleName': 'sample', 'RegistrationCode': 'person', 'Date': 'time', 'alloc': 'group'}

        # loading from Lab Data
        if type(data) == str:
            if data == 'microbiome':
                lab_data = MBLoader.MBLoader().get_data('segata_species', study_ids=self.params.study)
                fig_abundance_reads(lab_data.df_metadata)
            elif data == 'blood':
                lab_data = BloodTestsLoader.BloodTestsLoader().get_data(study_ids=self.params.study)
            elif data == 'body':
                lab_data = BodyMeasuresLoader.BodyMeasuresLoader().get_data(study_ids=self.params.study)
            elif data == 'cgm':
                lab_data = CGMLoader.CGMLoader().get_data(study_ids=self.params.study)
            else:
                raise Exception('data string is not valid')
            # TODO: check additional arguments for each loader

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
                file_df = file

            # unrecognizable file type
            else:
                raise Exception('file type is not valid')

            # fix known problems
            if 'RegistrationCode' in file_df.columns:
                file_df['RegistrationCode'] = file_df['RegistrationCode'].astype(str)

            if 'MeetingDate' in file_df.columns:
                file_df['MeetingDate'] = file_df['MeetingDate'] \
                    .dt.tz_localize('Asia/Jerusalem').dt.tz_convert('UTC')

            if 'Timestamp' in file_df.columns:
                file_df['Timestamp'] = pd.to_datetime(file_df['Timestamp'], errors='coerce') \
                    .dt.tz_localize('Asia/Jerusalem').dt.tz_convert('UTC')

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
            df1['log'] = False
            df2 = pd.DataFrame(np.log10(abundance_df).stack()).reset_index()
            df2['log'] = True
            df3 = pd.DataFrame(self.get_delta_df(abundance_df, self.params.controls['time_point']).stack())\
                .reset_index()
            df3['log'] = False
            df4 = pd.DataFrame(self.get_delta_df(np.log10(abundance_df), self.params.controls['time_point']).stack())\
                .reset_index()
            df4['log'] = True

            full_df = pd.concat([df1, df2, df3, df4]).rename(columns={0: 'abundance'})

            g = sns.FacetGrid(full_df, row='log', col='time_point', hue='group', palette=self.params.colors,
                              sharex=False, sharey=False, margin_titles=True)
            g = g.map(sns.distplot, 'abundance', hist=True, kde=False)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            plt.suptitle('Abundance distribution', horizontalalignment='left')
            plt.subplots_adjust(top=0.9)

            plt.savefig(os.path.join(self.dirs.figs, 'abundance distribution'))

        def fig_species_distribution(abundance_df):
            n_species_per_sample = pd.DataFrame(pd.DataFrame(
                (abundance_df > self.params.detection_threshold).sum(axis=1))
                                                .stack()).reset_index().rename(columns={0: 'n_species'})

            g = sns.FacetGrid(n_species_per_sample, col='time_point', hue='group', palette=self.params.colors,
                              sharex=True, sharey=True, margin_titles=True)
            g = g.map(sns.distplot, 'n_species', hist=True, kde=False)
            g.set_titles(col_template="{col_name}")
            for ax in g.axes.flatten():
                ax.set_title(label=ax.get_title(), color=self.params.colors[ax.get_title()])
            g.add_legend()
            g.set_ylabels('n_samples')
            plt.suptitle('Species distribution', horizontalalignment='left')
            plt.subplots_adjust(top=0.8)

            plt.savefig(os.path.join(self.dirs.figs, 'species distribution'))

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
            df = df.groupby(['group', 'person', 'time_point']).apply(lambda g: g.mean())
            # this obliviates any index that is not the group by!!!

        # filter out cases where a person has the wrong number of samples
        if n_entries is not None:
            df = df.groupby('person').filter(lambda x: len(x) == n_entries)

        # replace indices values according to dictionary
        if indices_dict is not None:
            for ind in indices_dict.keys():
                if ind in df.index.names:
                    df = df.rename(index=indices_dict[ind], level=ind)

        if obj.type == 'abundance':
            # fill missing values caused by MBLoader with detection threshold
            df = df.fillna(self.params.detection_threshold)

            # filter out empty species before declare_missing_values to shorten run time
            df = df.loc[:, (df != self.params.detection_threshold).any()]

            # declare missing values
            df = df.groupby(['person']).apply(declare_missing_values)

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
    def comp_stats(self, obj, test, between, delta=False, minimal_samples=0.1):
        """
        Compute the statistical significance for the difference between elements

        :param obj: (object) object to calculate the statistics on
        :param test: (str) mannwhitneyu, ttest_ind, ttest_rel, ttest_1samp, binom_test
        :param between: (str) group or time_point
        :param delta: (bool) whether or not to calculate the statistics based on the difference between
        time_point values
        :param minimal_samples: (percentage or int)

        :return: None
        """

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
                                     ncols=max(len(major_elements), len(minor_elements)))

                plt.subplots_adjust(bottom=0.5, left=0.14)
                plt.suptitle('{} - significant results'.format(test))

            ax = axes_internal[major_elements.index(major_e) * len(minor_elements) + minor_elements.index(minor_e)]

            # in case there are some significant results
            curr_data_df = curr_data_df.dropna()
            if curr_data_df.shape[0] != 0:

                # transfer the columns to rows with identifiers
                curr_data_df = pd.DataFrame(curr_data_df.stack()).reset_index() \
                    .rename(columns={'level_0': 'index', 0: 'list'})

                # split the lists to a row for each element in them
                curr_data_df = curr_data_df \
                    .merge(curr_data_df['list'].apply(pd.Series), right_index=True, left_index=True) \
                    .melt(id_vars=['index', major, 'list'], value_name='value') \
                    .drop(['list', 'variable'], axis=1) \
                    .dropna()

                # clean the bacteria name
                if obj.columns == 'bacteria':
                    curr_data_df['index'] = curr_data_df['index'].apply(
                            lambda row: '{} ({})'.format(row.split('|')[6], row.split('|')[-1]).replace('s__', ''))
                    y_label = 'abundance'
                else:
                    y_label = '{} values'.format(obj.columns)

                # clean the blood test name
                if obj.type == 'blood':
                    curr_data_df['index'] = curr_data_df['index'].apply(
                        lambda row: row.replace('bt__', ''))

                # sort index values so they will be sorted in the plot
                curr_data_df = curr_data_df.sort_values(by='index')

                # plotting
                sns.boxplot(x='index', y='value', hue=major, palette=self.params.colors, data=curr_data_df, ax=ax)
                ax.set_ylabel(y_label)
                ax.set_xlabel('')
                ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=90)
                ax.legend().set_visible(False)

            else:
                # texting
                ax.text(x=0.5, y=0.5, s='no significant results', horizontalalignment='center')
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
                    'stats {}{} {} - {}'.format(obj.type, delta_str, test, stats_df_list[-1].name.replace('.', ''))))
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
        major_elements = df.index.get_level_values(major).unique().values.tolist()
        minor_elements = df.index.get_level_values(minor).unique().values.tolist()

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
                if test not in tests_without_control:
                    stats_df_list[-1].name = '{} {} vs. {} of {} {}'.format(major, major_control, major_e,
                                                                            minor, minor_e)
                    excel_sheet = '{} {} {}'.format(major_control, major_e, minor_e).replace('mediterranean', 'med')
                    # mediterranean handling just to limit excel sheet name to the 31 character count limit
                else:
                    stats_df_list[-1].name = '{} {} of {} {}'.format(major, major_e, minor, minor_e)
                    excel_sheet = '{} {}'.format(major_e, minor_e)

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
                    if test == 'ttest_rel':
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

                        # this is only relevant when between is 'time'

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

        # excel finishes
        excel_writer.save()
        excel_writer.close()

    def score_models(self, xobj, yobj, delta=True, minimal_samples=0.9,
                     model_type='linear', n_iterations=10, n_splits=5, random_seed=None):
        """
        Calculate n_iterations prediction models for the yobj based on the xobj
        and return the summary statistics on the models scores

        :param xobj: (object) the features object
        :param yobj: (object) the target object
        :param delta: (bool) whether or not to calculate the models based on the difference between time_point values
        :param minimal_samples: (percentage or int) that needs to be included in the models
        and therefore cannot be deleted in the missing values filtering process
        :param model_type: (str) 'linear' or 'xgb'
        :param n_iterations: (int) number of models to be created for each target
        :param n_splits: (int) number of time to split the data in the KFold
        :param random_seed: (int) seed for the random process used for splitting the data

        :return: (pd.DataFrame) summary_score_df
        """

        # retrieving the data frames from the objects
        if type(xobj) == _Object and type(yobj) == _Object:
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
        if xobj.type == 'abundance':  # only for abundance otherwise the min min might not be the best solution
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
        x = np.array(x_df)

        # create a data frame to fill with results
        params = ['R^2', 'r', 'p', 'roc_auc']
        columns = pd.MultiIndex.from_product([np.arange(n_iterations), params], names=['iteration', 'statistic'])
        score_df = pd.DataFrame(index=y_df.columns, columns=columns)

        # random_state is used by the xgb model to split the data the same way in the KFold
        if random_seed is not None:
            np.random.seed(random_seed)
        random_state = np.random.randint(100, size=n_iterations)  # the 100 is arbitrary

        # create a model for each target column
        for col in y_df.columns:

            # check if column is binary in order to know which model to use
            if sorted(y_df[col].unique()) == [0, 1] or sorted(y_df[col].unique()) == [False, True]:
                binary = True
            else:
                binary = False

            # skip column in some cases
            if model_type == 'linear' and binary:
                print('skipping binary column while model is linear - {}'.format(col))
                continue
            if y_df[col].isna().sum() != 0:
                print('skipping column with {} missing values - {}'.format(y_df[col].isna().sum(), col))
                continue
            if y_df[col].unique().shape[0] == 1:
                print('skipping column with a single value ({}) - {}'.format(y_df[col].unique(), col))
                continue

            # conversion
            y = np.array(y_df[col])

            # run iterations of the model for scoring statistics
            for i in np.arange(n_iterations):

                # empty array to hold all predictions in order to score them later together
                y_pred = np.empty(y.shape)
                y_pred.fill(np.nan)

                # create K folds of the data and do for each fold
                kf = KFold(n_splits=n_splits, random_state=random_state[i], shuffle=True)
                for train_index, test_index in kf.split(x):

                    # define the subsets
                    x_train, x_test = x[train_index], x[test_index]
                    y_train, y_test = y[train_index], y[test_index]

                    # linear model
                    if model_type == 'linear':
                        model = LinearRegression().fit(x_train, y_train)

                    # xgb model
                    elif model_type == 'xgb' and binary:
                        model = XGBClassifier().fit(x_train, y_train)
                    elif model_type == 'xgb' and not binary:
                        model = XGBRegressor(objective='reg:squarederror').fit(x_train, y_train)
                        # objective='reg:squarederror' is the default and is written explicitly just to avoid a warning
                    # TODO: ask Eran about the parameters
                    else:
                        raise Exception('model not valid')

                    # fill up y_pred in steps so can be analyzed in total later on
                    y_pred[test_index] = model.predict(x_test)

                # produce model evaluation
                if binary:
                    score_df.loc[col, (i, 'roc_auc')] = roc_auc_score(y, y_pred)
                else:
                    score_df.loc[col, (i, 'R^2')] = r2_score(y, y_pred)
                    # R^2 (coefficient of determination) regression score function.
                    # Best possible score is 1.0 and it can be negative (because the model can be arbitrarily worse).
                    # A constant model that always predicts the expected value of y, disregarding the input features,
                    # would get a R^2 score of 0.0.
                    score_df.loc[col, (i, ['r', 'p'])] = pearsonr(y, y_pred)
                    # The p-values are not entirely reliable but are probably reasonable for data sets larger than 500

        # dropping all the un-ran columns
        score_df = score_df.dropna(axis=0, how='all')

        # summarize the scores using the iterations scores distribution
        columns = pd.MultiIndex.from_product([params, ['mean', 'std']], names=['statistic', 'distribution'])
        summary_score_df = pd.DataFrame(index=score_df.index, columns=columns)

        for param in params:
            summary_score_df.loc[:, (param, 'mean')] = score_df.xs(param, level='statistic', axis=1).mean(axis=1)
            summary_score_df.loc[:, (param, 'std')] = score_df.xs(param, level='statistic', axis=1).std(axis=1)

        summary_score_df = summary_score_df.sort_values((params[0], 'mean'), ascending=False)

        # save the data frame as excel
        # create excel writer to fill up
        excel_path = os.path.join(self.dirs.excels,
                                  'models {} {}{} {}.xlsx'.format(xobj.type, yobj.type, delta_str, model_type))
        excel_writer = pd.ExcelWriter(excel_path)
        summary_score_df.to_excel(excel_writer, freeze_panes=(2, 1))
        excel_writer.save()
        excel_writer.close()

        return summary_score_df

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

        diversity_df = pd.DataFrame(-(abundance_df * np.log10(abundance_df)).sum(axis=1).dropna())
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

    # general file handling
    @staticmethod
    def ftp_download(address, username, password,
                     directories, skip_files=None, destination=os.getcwd(),
                     check_only=True):
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

    @staticmethod
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


# Helper classes
class _Object:

    # Initialization
    def __init__(self, obj_type=None, df=None, columns=None):
        self.type = obj_type
        self.df = df
        self.columns = columns


class _Objects:

    # Initialization
    def __init__(self):

        # microbiome
        self.abundance = _Object(obj_type='abundance', columns='bacteria')
        self.diversity = _Object(obj_type='diversity', columns='diversity')

        # others
        self.blood = _Object(obj_type='blood', columns='tests')
        self.body = _Object(obj_type='body', columns='measurements')


class _Directories:

    # Initialization
    def __init__(self, base_directory):

        self.figs = os.path.join(base_directory, 'figs')
        self.excels = os.path.join(base_directory, 'excels')
        self.data_frames = os.path.join(base_directory, 'data_frames')

        directories = [self.figs, self.excels, self.data_frames]

        for Dir in directories:
            if not os.path.exists(Dir):
                os.mkdir(Dir)


class _Parameters:

    # Initialization
    def __init__(self, study=None, controls=None, colors=None, alpha=None, detection_threshold=None):

        self.study = study
        self.controls = controls
        self.colors = colors
        self.alpha = alpha
        self.detection_threshold = detection_threshold


if __name__ == "__main__":

    print(help(Study))
