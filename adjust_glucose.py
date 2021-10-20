# import time
# import os
# from Utils import Load, Write, is_file_newer, get_file_dates, isnumeric
# import Defs
# # from HighPPGRs import HighPPGRs
# import pandas as pd
# from pandas import Series, DataFrame, concat
# import numpy as np
# from itertools import product
# import datetime
#
# class AnalysisData(object):
#     def __init__(self, user = ''):
#      #   self._external_file_dates()
#
#         self._randomization = Defs._clear_and_get_cached(Defs._united_randomization_cache, self._get_united_randomization_file, cache_loc=Defs._cache_dir)
#         self._meal_ppgrs = Defs._clear_and_get_cached(Defs._meal_ppgrs_cache, Defs._read_excel, (Defs._meal_ppgrs_file,), cache_loc=Defs._cache_dir)
#         self._blood_tests = Defs._clear_and_get_cached(Defs._blood_tests_cache, self._load_blood_tests, cache_loc=Defs._cache_dir)
#
#         self._cgm = Defs._clear_and_get_cached(Defs._cgm_cache, self._load_cgm_data, cache_loc=Defs._cache_dir)
#         print 'Connections={}'.format(len(self._cgm.index.get_level_values('ConnectionID').unique()))
#
#         self._adjusted_glucose = Defs._clear_and_get_cached(Defs._adjusted_glucose_cache, self._adjust_cgm, cache_loc=Defs._cache_dir)
#
#         self._cgm_inactive = Defs._clear_and_get_cached(Defs._cgm_inactive_cache, self._load_cgm_data, (0,), cache_loc=Defs._cache_dir)
#         self._adjusted_glucose_inactive = Defs._clear_and_get_cached(Defs._adjusted_glucose_inactive_cache, self._adjust_cgm, (0,), cache_loc=Defs._cache_dir)
#
#         self._con_metadata, self._user_metadata = Defs._clear_and_get_cached(Defs._user_metadata_cache, self._gen_user_metadata, cache_loc=Defs._cache_dir)
#         self._glucose_features = Defs._clear_and_get_cached(Defs._glucose_features_cache, self._gen_glucose_thr_features, cache_loc=Defs._cache_dir)
#         self._meals_df = Defs._clear_and_get_cached(Defs._meals_df_cache, self._gen_meals_df, cache_loc=Defs._cache_dir)
#         self._dietary_features = Defs._clear_and_get_cached(Defs._dietary_features_cache, self._gen_dietary_features, cache_loc=Defs._cache_dir)
#         self._users_plot_dfs = Defs._clear_and_get_cached(Defs._prepared_user_plots_cache, self._prepare_user_plot_dfs, cache_loc=Defs._cache_dir)
#         self._outcome_stats = Defs._clear_and_get_cached(Defs._outcome_stats_cache, self._get_outcome_stats, cache_loc=Defs._cache_dir)
#         self.user = user
#         # self._extract_all_hba1c_dates()
#
#     def _filter_study_df(self, df, k, user_ids, reg_ids):
#         org_size = len(df)
#         res = df
#         if df.index.names[0] == 'UserID':
#             res = df[df.index.get_level_values(0).isin(user_ids)]
#         elif df.index.names[0] == 'RegistrationCode':
#             res = df[df.index.get_level_values(0).isin(reg_ids)]
#         elif 'UserID' in df.columns:
#             res = df[df['UserID'].isin(user_ids)]
#         elif 'RegistrationCode' in df.columns:
#             res = df[df['RegistrationCode'].isin(reg_ids)]
#         else:
#             print 'No keys for {} type {}'.format(k, type(df))
#         print 'Reduced {} from {} to {}'.format(k, org_size, len(res))
#         return res
#
#     def _filter_study(self, study):
#         print 'Selecting study {}'.format(study)
#         reg_ids = self._randomization[self._randomization['Study'] == study]['RegistrationCode'].values
#         user_ids = np.array(list(self._randomization[self._randomization['Study'] == study].index))
#         for k in self.__dict__.keys():
#             df = self.__dict__[k]
#             if type(df) == DataFrame:
#                 self.__dict__[k] = self._filter_study_df(df, k, user_ids, reg_ids)
#             elif k == '_users_plot_dfs':
#                 for k in df.__dict__.keys():
#                     subdf = df.__dict__[k]
#                     if type(subdf) == DataFrame:
#                         df.__dict__[k] = self._filter_study_df(subdf, k, user_ids, reg_ids)
#                     else:
#                         print 'Skipping filtering of {}'.format(k)
#             else:
#                 print 'Skipping filtering of {}'.format(k)
#
#     def _batch_num(self, reg_id):
#         return int(self._user_metadata[self._user_metadata['RegistrationCode'] == reg_id]['BatchNum'].values[0])
#
#     def _is_main(self, reg_id):
#         return self._user_metadata[self._user_metadata['RegistrationCode'] == reg_id]['is main'].values[0] == 1
#
#     def _user_id_2_reg_id(self, user_id):
#         return self._user_metadata.loc[user_id]['RegistrationCode']
#
#     def _reg_id_2_user_id(self, reg_id):
#         return self._user_metadata[self._user_metadata['RegistrationCode']==reg_id].index.get_level_values(0)[0]
#
#     def _reg_id_2_gender(self, reg_id):
#         return int(self._user_metadata[self._user_metadata['RegistrationCode']==reg_id]['Gender'].values[0])
#
#     def _get_user_study(self, reg_id):
#         return self._randomization.loc[self._reg_id_2_user_id(reg_id)]['Study']
#
#     def _get_user_ids_for_study(self, study):
#         return list(self._randomization[self._randomization['Study']==study].index.get_level_values('UserID'))
#
#     def _get_reg_ids_for_study(self, study, only_active=False):
#         if only_active:
#             return self._randomization[np.logical_and(self._randomization['Study'] == study, self._randomization['IsActive'])]['RegistrationCode'].values
#         else:
#             return self._randomization[self._randomization['Study']==study]['RegistrationCode'].values
#
#     def _reg_ids_by_days_from_t0(self, min_days_from_t0=-np.inf, max_days_from_t0=np.inf):
#         return self._user_metadata[np.logical_and(self._user_metadata['DaysFromT0']>=min_days_from_t0,
#                                                   self._user_metadata['DaysFromT0']<=max_days_from_t0)]
#
#     def _filter_by_feature_by_start_month(self, df, feature, start_month):
#         assert feature in Defs._feature_2_min_days
#         if type(df)==Series:
#             df = df.to_frame(df.name)
#         df = df.join((pd.Timestamp(datetime.date.today())-self._randomization['T0']).apply(lambda x: x.days).rename('DaysFromT0'))
#         # df = df.join(self._user_metadata['DaysFromT0'])
#         return df[df['DaysFromT0']>=(start_month-1)*30+Defs._feature_2_min_days[feature]].iloc[:,:-1]
#
#     def _get_outcome_stats(self):
#         p = self._users_plot_dfs
#
#         m_f = lambda x: x.iloc[:,1]-x.iloc[:,0]
#
#         def _add_glucose_feature(glucose_adj_percentile_name, feature, group, threshold, value, months):
#             prefix = [glucose_adj_percentile_name, feature, group, threshold]
#             base_str = '{} {} {} {}'.format(glucose_adj_percentile_name, feature, group, threshold)
#             val_str = '{} {} {} {} {} vs Prf'.format(glucose_adj_percentile_name, feature, group, threshold, value)
#             cmp_f = Defs._time_in_range_cmp_p if feature.endswith('_p') else Defs._time_in_range_cmp
#             f = lambda x: cmp_f(self._filter_by_feature_by_start_month(x.dropna(), feature, value))
#             return [(tuple(prefix + [0]), None, base_str + ' 0M', 'Glucose', None),
#                     (tuple(prefix + [value]), None, base_str + ' {}M'.format(months), 'Glucose', None),
#                     (tuple(prefix + [0]), tuple(prefix + [value]), val_str, 'Glucose', f)]
#
#         def _add_blood_test_feature(feature):
#             return [((feature, 'T0'), None, '{} 0M'.format(feature), 'Blood', None),
#                     ((feature, '3M'), None, '{} 3M'.format(feature), 'Blood', None),
#                     ((feature, '6M'), None, '{} 6M'.format(feature), 'Blood', None),
#                     ((feature, 'T0'), (feature, '3M'), '{} 3M vs T0'.format(feature), 'Blood', m_f),
#                     ((feature, 'T0'), (feature, '6M'), '{} 6M vs T0'.format(feature), 'Blood', m_f)]
#
#         def _add_max_blood_test_feature(feature):
#             return [((feature, 'Max_Screening_T0'), None, '{} 0M'.format(feature), 'Blood', None),
#                     ((feature, '3M'), None, '{} 3M'.format(feature), 'Blood', None),
#                     ((feature, '6M'), None, '{} 6M'.format(feature), 'Blood', None),
#                     ((feature, 'Max_Screening_T0'), (feature, '3M'), '{} 3M vs MaxScrT0'.format(feature), 'Blood', m_f),
#                     ((feature, 'Max_Screening_T0'), (feature, '6M'), '{} 6M vs MaxScrT0'.format(feature), 'Blood', m_f)]
#
#         def _add_max_weight_feature(feature):
#             return [((feature, 'Max_Profiling_T0'), None, '{} 0M'.format(feature), 'Weight', None),
#                     ((feature, '3M'), None, '{} 3M'.format(feature), 'Weight', None),
#                     ((feature, '6M'), None, '{} 6M'.format(feature), 'Weight', None),
#                     ((feature, 'Max_Profiling_T0'), (feature, '3M'), '{} 3M vs MaxPrfT0'.format(feature), 'Weight', m_f),
#                     ((feature, 'Max_Profiling_T0'), (feature, '6M'), '{} 6M vs MaxPrfT0'.format(feature), 'Weight', m_f)]
#
#         def _add_diet_3M_feature(feature):
#             return [((feature, 0), None, '{} 0M'.format(feature), 'Diet_3M', None),
#                     ((feature, 1), None, '{} 3M'.format(feature), 'Diet_3M', None),
#                     ((feature, 0), (feature, 1), '{} 3M vs T0'.format(feature), 'Diet_3M', m_f)]
#
#         def _add_diet_6M_feature(feature):
#             return [((feature, 1), None, '{} 6M'.format(feature), 'Diet_6M', None),
#                     ((feature, 0), (feature, 1), '{} 6M vs T0'.format(feature), 'Diet_6M', m_f)]
#
#         df_name2df = {'Blood': p.blood_tests_table,
#                       'Weight': p.weights_df_table,
#                       'Glucose': self._glucose_features,
#                       'Diet_3M': self._dietary_features.by_3months.unstack('3MonthsFromT0'),
#                       'Diet_6M': self._dietary_features.by_6months.unstack('6MonthsFromT0'),}
#
#         outcomes = []
#         for feature in ['Carbs % kcal', 'Sodium', 'TransFat % kcal', 'MonoUnsaturatedFat % kcal', 'SaturatedFat % kcal',
#                         'Fiber g/1000 kcal', 'Lipid % kcal', 'SugarsTotal g/1000 kcal', 'PolyUnsaturatedFat % kcal', 'Protein % kcal',
#                         'Kcals per day', 'Kcals per day % of target', 'main_score', 'control_score_carbs', 'control_score_fat',
#                         'control_score_saturated_fat', 'control_score_fiber', 'control_score_macronutrients', 'control_score']:
#             outcomes += _add_diet_3M_feature(feature)
#             outcomes += _add_diet_6M_feature(feature)
#         for feature in ['HbA1C%', 'FPG']:
#             outcomes += _add_max_blood_test_feature(feature)
#         for feature in ['Insulin', 'Triglycerides', 'Fructosamine', 'Cholesterol - HDL', 'LDL Cholesterol']:
#             outcomes += _add_blood_test_feature(feature)
#         for feature in ['Weight', 'Fat_perc', 'BP_sys', 'BP_dia', 'BMI', 'BMR', 'HeartRate', 'Waist', 'Hips']:
#             outcomes += _add_max_weight_feature(feature)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By3Months', 'ByConnection', '130', 1, 3)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By3Months', 'ByConnection', '140', 1, 3)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By6Months', 'ByConnection', '130', 1, 6)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By6Months', 'ByConnection', '140', 1, 6)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By3Months', 'ByDay', 'Mean', 1, 3)
#         outcomes += _add_glucose_feature(Defs.GLUCOSE_ADJ_COL, 'By6Months', 'ByDay', 'Mean', 1, 6)
#         outcomes += _add_glucose_feature('PPGR', 'By3Months', 'ByConnection', '95_p', 1, 3)
#         outcomes += _add_glucose_feature('PPGR', 'By6Months', 'ByConnection', '95_p', 1, 6)
#
#         res = []
#         for base_col, out_col, out_name, df_name, f in outcomes:
#             if out_col is None:
#                 out = df_name2df[df_name][base_col]
#             else:
#                 df = concat([df_name2df[df_name][base_col], df_name2df[df_name][out_col]], 1)
#                 out = f(df)
#             if df_name == 'Glucose':
#                 out = Defs._change_index(out.to_frame(out_name), self._user_metadata, 'RegistrationCode').sort_index()
#             else:
#                 out.name = out_name
#             res.append(out)
#
#         res = concat(res, 1)
#         res = res.join(self._user_metadata[['RegistrationCode', 'is main', 'Study', 'Calories']].set_index('RegistrationCode')).\
#             join(self._randomization.set_index('RegistrationCode')[['Age', 'Gender']]).\
#             sort_values('HbA1C% 3M vs MaxScrT0')
#         res['Colors'] = 'red'
#         res.loc[res[res['is main'] == 1].index, 'Colors'] = 'green'
#         # res = res[~np.isnan(res['HbA1C% 3M vs MaxScrT0'])]
#         return res
#
#     def _get_baseline_stats(self):
#         res = []
#         cgm_thresholds = [(t, False) for t in Defs._cgm_high_thresholds] + [(t, True) for t in Defs._cgm_high_p_thresholds]
#         for cgm_threshold, is_dynamic in cgm_thresholds:
#             for glucose_adj_percentile_name in ['PPGR', 'GlucoseAdj50']:
#                 if glucose_adj_percentile_name == 'PPGR' and not is_dynamic:
#                     continue
#                 features = self._adjusted_glucose[[glucose_adj_percentile_name, 'DaysFromT0']]
#                 features = features[features['DaysFromT0']<0][glucose_adj_percentile_name].dropna().groupby(level=['UserID'])
#                 features = features.apply(lambda x: np.percentile(x, cgm_threshold)) if is_dynamic \
#                     else features.apply(lambda x: Defs._time_in_range(x, cgm_threshold))
#                 features.name = '{}_{}{}'.format(glucose_adj_percentile_name, cgm_threshold, '_p' if is_dynamic else '')
#                 res.append(features)
#         res = concat(res, 1)
#         res = res.join(self._user_metadata[['DaysFromT0', 'NumProfilingDays']])
#         res = Defs._change_index(res, self._user_metadata, 'RegistrationCode').sort_index()
#
#         profiling_blood_tests = self._users_plot_dfs.blood_tests_df[self._users_plot_dfs.blood_tests_df.MeetingTypeID.apply(lambda x: x in ['T0', 'Screening'])]
#         profiling_blood_tests = profiling_blood_tests.set_index(['RegistrationCode', 'BloodTest', 'MeetingTypeID']).iloc[:,0].apply(pd.to_numeric)
#         profiling_blood_tests = profiling_blood_tests.groupby(['RegistrationCode', 'BloodTest', 'MeetingTypeID']).mean().unstack(
#             ['BloodTest', 'MeetingTypeID'])
#         self._fix_blood_tests(profiling_blood_tests)
#
#         valid_weights_df = self._users_plot_dfs.weights_df[np.logical_and(self._users_plot_dfs.weights_df['Weight']>0, self._users_plot_dfs.weights_df['Fat_perc']>0)]
#         profiling_weights = valid_weights_df[valid_weights_df.MeetingTypeID.apply(lambda x: x in ['T0', 'Profiling'])] \
#             [['RegistrationCode', 'MeetingTypeID', 'Weight', 'Fat_perc']].set_index(['RegistrationCode', 'MeetingTypeID']).\
#             unstack('MeetingTypeID')
#         for weight_feature in ['Weight']:
#             profiling_weights[(weight_feature, 'Max_Profiling_T0')] = profiling_weights[[(weight_feature, 'T0'), (weight_feature, 'Profiling')]].max(1)
#             profiling_weights[(weight_feature, 'Avg_Profiling_T0')] = profiling_weights[[(weight_feature, 'T0'), (weight_feature, 'Profiling')]].mean(1)
#
#         res = concat([res, profiling_blood_tests, profiling_weights], 1)
#         return res
#
#     def _get_compliance_stats(self):
#         org_df = self._users_plot_dfs.score_stats_df.set_index(['RegistrationCode', 'MonthsFromT0'])
#         carbs_perc = self._users_plot_dfs.by_month_carbs_fat_cals_by_meals_df[self._users_plot_dfs.by_month_carbs_fat_cals_by_meals_df['Meal score'] == 'Carb%'].set_index('RegistrationCode')\
#             [['MonthsFromT0', 'Energy']].rename({'Energy': 'Carbs%'}, axis=1).set_index('MonthsFromT0', append=True)
#         energy_by_time = self._users_plot_dfs.by_month_energy_by_time.set_index(['RegistrationCode', 'MonthsFromT0'])
#         org_df = org_df.join(carbs_perc).join(energy_by_time).reset_index('MonthsFromT0')
#
#         res = [org_df]
#         avg_months = [(1,3), (4,6), (1,6)]
#         for (month_start, month_end) in avg_months:
#             df = org_df.loc[list(org_df[org_df.MonthsFromT0==month_end].index.get_level_values(0))]
#             df = df[np.logical_and(df.MonthsFromT0>=month_start, df.MonthsFromT0<=month_end)].drop('MonthsFromT0', 1).groupby(level=['RegistrationCode']).mean()
#             df['MonthsFromT0'] = '{}M-{}M'.format(month_start, month_end)
#             res.append(df)
#
#         res = concat(res).set_index('MonthsFromT0', append=True).unstack('MonthsFromT0')
#         res.columns.names = ['ScoreFeature', 'MonthsFromT0']
#
#         return res
#
#     def _gen_user_metadata(self):
#         active_users = Series(index=self._randomization[self._randomization['IsActive'] == 1].index, name='ActiveUsers').apply(lambda x: 1)
#         con_header = self._cgm.groupby(['UserID', 'ConnectionID']).head(1).reset_index('GlucoseTimestamp', drop=True).drop('T0', axis='columns').join(self._randomization['T0'])
#         is_profiling = con_header['Period_end'] < con_header['T0']
#         con_end_day_after_T0 = (con_header['Period_end'] - con_header['T0']).apply(lambda x: x.days)
#         con_start_day_after_T0 = (con_header['Period_start'] - con_header['T0']).apply(lambda x: x.days)
#         con_mid_day_after_T0 = (con_end_day_after_T0 + con_start_day_after_T0) / 2
#         con_metadata = self._cgm[['RegistrationCode', 'is main']]. \
#             reset_index('GlucoseTimestamp', drop=True).reset_index(['UserID', 'ConnectionID']). \
#             drop_duplicates(['UserID', 'ConnectionID']).set_index(['UserID', 'ConnectionID']). \
#             join(is_profiling.to_frame('is profiling')). \
#             join(self._randomization[['T0', 'BatchNum', 'Study', 'Gender']]). \
#             join(con_start_day_after_T0.to_frame('StartDayAfterT0')). \
#             join(con_end_day_after_T0.to_frame('EndDayAfterT0')). \
#             join(con_mid_day_after_T0.to_frame('MidDayAfterT0')). \
#             join(con_header[['Period_start', 'Period_end']]). \
#             join(con_end_day_after_T0.apply(Defs._days_after_T0_2_meeting_type).to_frame('MeetingType'))
#
#         cons_per_user = self._adjusted_glucose.reset_index('ConnectionID').drop_duplicates('ConnectionID').groupby(level='UserID').size()
#         cons_per_user.name = 'ConsPerUser'
#         users_3M = self._blood_tests[self._blood_tests['MeetingTypeID'] == 9]['UserID'].drop_duplicates()
#         users_3M = Series(np.ones(users_3M.shape[0]), index=users_3M, name='3M')
#         users = con_metadata[['RegistrationCode', 'T0', 'BatchNum', 'Study', 'is main', 'Gender']].reset_index('ConnectionID', drop=True).drop_duplicates()
#         print cons_per_user.shape
#         print active_users.shape
#         print users.shape
#         print  users_3M.shape
#         user_metadata = concat([cons_per_user, users_3M, active_users, users], 1)
#         user_metadata = user_metadata.loc[user_metadata.reset_index()['UserID'].dropna()]
#         user_metadata = user_metadata.join(self._adjusted_glucose.groupby('UserID').apply(lambda x: len(x[x['DaysFromT0'] < 0]['DaysFromT0'].unique())).to_frame('Profiling Days'))
#         user_metadata = user_metadata.join(self._randomization[['Dietitian', 'Result']])
#         user_metadata['DaysFromT0'] = (pd.Timestamp(datetime.date.today()) - user_metadata['T0']).apply(lambda x: x.days)
#         user_metadata['NumProfilingDays'] = self._adjusted_glucose[self._adjusted_glucose['DaysFromT0'] < 0].groupby(level='UserID').size() / 96
#         user_metadata = user_metadata.join(con_metadata[con_metadata['is profiling']].sort_values('EndDayAfterT0', ascending=False).groupby(level=['UserID']).first()
#                                            [['EndDayAfterT0', 'Period_end']].\
#                                            rename({'EndDayAfterT0': 'ProfilingEndByDayFromT0', 'Period_end': 'ProfilingEnd'}, axis='columns'))
#         user_metadata['ProfilingStartByDayFromT0'] = con_metadata.groupby('UserID')['StartDayAfterT0'].min().sort_values()
#         user_metadata.rename({'Result': 'Calories'}, axis=1, inplace=True)
#         return con_metadata, user_metadata
#
#     def _gen_glucose_mean_feature(self, df, groupby_f, groupby_num_days):
#         return df.groupby(df['DaysFromT0'].apply(groupby_f)).apply(\
#             lambda x: Defs._apply_with_len_filter(x, lambda x: x.mean(), groupby_num_days)).\
#             drop('DaysFromT0', 1).dropna()
#
#     def _gen_glucose_feature(self, df, groupby_f, groupby_num_days, cgm_threshold, is_dynamic_threshold):
#         if is_dynamic_threshold:
#             if len(df[df['DaysFromT0'] < 0].dropna()) == 0:
#                 return df.groupby(df['DaysFromT0'].apply(groupby_f)).apply(lambda x: x.apply(lambda x: np.NaN)).drop('DaysFromT0', 1).dropna()
#             cgm_threshold = np.percentile(df[df['DaysFromT0']<0].iloc[:,0].dropna(), cgm_threshold)
#         return df.groupby(df['DaysFromT0'].apply(groupby_f)).apply(\
#             lambda x: Defs._apply_with_len_filter(x, lambda x: Defs._time_in_range(x, cgm_threshold), groupby_num_days)).\
#             drop('DaysFromT0', 1).dropna()
#
#     def _gen_glucose_feature_by_con(self, df, groupby_f, groupby_num_days, cgm_threshold, is_dynamic_threshold):
#         if is_dynamic_threshold:
#             if len(df[df['MidDayAfterT0'] < 0].dropna()) == 0:
#                 return df.groupby(df['MidDayAfterT0'].apply(groupby_f)).apply(lambda x: x.apply(lambda x: np.NaN)).drop('MidDayAfterT0', 1).dropna()
#             cgm_threshold = df[df['MidDayAfterT0']<0].iloc[:,0].dropna().groupby(level='ConnectionID').apply(lambda x: np.percentile(x, cgm_threshold)).median()
#         con_group = df['MidDayAfterT0'].groupby(level='ConnectionID').head(1).apply(groupby_f).reset_index(['UserID', 'GlucoseTimestamp'], drop=True)
#         sufficient_cgm_vals = df.iloc[:, 0].groupby(level='ConnectionID').count().to_frame('Count').join(con_group).groupby('MidDayAfterT0').sum().apply(lambda x: Defs._sufficient_cgm_values(x, groupby_num_days))
#         res = df.iloc[:,0].groupby(level='ConnectionID').apply(lambda x: Defs._time_in_range(x, cgm_threshold)).to_frame(df.columns[0]).join(con_group).groupby('MidDayAfterT0').median().dropna().join(sufficient_cgm_vals)
#         return res[res['Count']].iloc[:,0].to_frame(res.columns[0])
#
#     def _gen_glucose_thr_features(self):
#         def _fix_features(features, cgm_threshold, is_dynamic, groupby_name, connection_grouping_name, name):
#             features['Threshold'] = '{}{}'.format(cgm_threshold, '_p' if is_dynamic else '')
#             features['Feature'] = groupby_name
#             features['Groupby'] = connection_grouping_name
#             features = features.set_index(['Feature', 'Groupby', 'Threshold'], append=True).unstack('Feature').unstack('Groupby').unstack('Threshold').unstack(name)
#             features.columns.names = ['GlucoseAdj', 'Feature', 'Groupby', 'Threshold', 'Value']
#             return features
#
#         gv = self._adjusted_glucose.reset_index('GlucoseTimestamp').join(self._con_metadata['MidDayAfterT0']).set_index('GlucoseTimestamp', append=True)
#         # groups = [(Defs._month_groups, 30, 'ByMonth'), (Defs._3month_groups, 90, 'By3Months'), (Defs._6month_groups, 180, 'By6Months'), (Defs._months2_3_groups, 60, 'ByMonths2_3'), (Defs._months2_6_groups, 150, 'ByMonths2_6')]
#         groups = [(Defs._month_groups, 30, 'ByMonth'), (Defs._3month_groups, 90, 'By3Months'), (Defs._6month_groups, 180, 'By6Months')]
#         cgm_thresholds = [(t, False) for t in Defs._cgm_high_thresholds] + [(t, True) for t in Defs._cgm_high_p_thresholds]
#         res = []
#         for groupby_f, groupby_num_days, groupby_name in groups:
#             for glucose_adj_percentile_name in Defs._glucose_adj_percentiles_neighbor_names + ['GlucoseValue']:
#                 print 'Generating mean feature {} for {}'.format(groupby_name, glucose_adj_percentile_name)
#                 features = gv[[glucose_adj_percentile_name, 'DaysFromT0']].dropna().groupby(level=['UserID']).apply(\
#                     lambda x: self._gen_glucose_mean_feature(x, groupby_f, groupby_num_days))
#                 res.append(_fix_features(features, 'Mean', False, groupby_name, 'ByDay', 'DaysFromT0'))
#         for cgm_threshold, is_dynamic in cgm_thresholds:
#             for groupby_f, groupby_num_days, groupby_name in groups:
#                 for glucose_adj_percentile_name in ['PPGR'] + Defs._glucose_adj_percentiles_names + Defs._glucose_adj_percentiles_neighbor_names:
#                     if glucose_adj_percentile_name == 'PPGR' and not is_dynamic:
#                         continue
#                     print 'Generating feature {} {} for {}'.format(groupby_name, cgm_threshold, glucose_adj_percentile_name)
#                     features = gv[[glucose_adj_percentile_name, 'MidDayAfterT0']].dropna().groupby(level=['UserID']).apply(\
#                         lambda x: self._gen_glucose_feature_by_con(x, groupby_f, groupby_num_days, cgm_threshold, is_dynamic))
#                     res.append(_fix_features(features, cgm_threshold, is_dynamic, groupby_name, 'ByConnection', 'MidDayAfterT0'))
#                     features = gv[[glucose_adj_percentile_name, 'DaysFromT0']].dropna().groupby(level=['UserID']).apply(\
#                         lambda x: self._gen_glucose_feature(x, groupby_f, groupby_num_days, cgm_threshold, is_dynamic))
#                     res.append(_fix_features(features, cgm_threshold, is_dynamic, groupby_name, 'ByDay', 'DaysFromT0'))
#         return concat(res, 1)
#
#     def _adjust_connection_cgm_by_percentile(self, x, percentile, p):
#         shift = np.percentile(x, p) - percentile
#         # print 'Fit Con={} p={} Percentile={} Median={} ShiftBy={:.3f}'.\
#         #     format(x.index.get_level_values('ConnectionID')[0], p, np.percentile(x, p), x.median(), shift)
#         return Series(x-shift, index=x.index)
#
#     def _adjust_connection_cgm_by_neighbors(self, x, p, n):
#         y_percentile = x.groupby('ConnectionID').apply(lambda x: np.percentile(x, p))
#         return x + y_percentile.rolling(n, min_periods=1, center=True).median() - y_percentile
#
#     def _adjust_user_cgm(self, x):
#         user_id = x.index.get_level_values('UserID')[0]
#         print 'Adjusting glucose values for UserID {}'.format(user_id)
#
#         gvs, gv_names = [], []
#
#         cons = self._get_user_connection_times_from_t0(user_id).reset_index('UserID', drop=True)
#         cons_intervention = cons[cons<250]
#         cons_post_intervention = cons[cons>=250]
#
#         def rolling(y_modes, n, cons_intervention, cons_post_intervention, is_median):
#             y_modes_intervention = y_modes.loc[cons_intervention.index.get_values()].rolling(n, min_periods=1, center=True)
#             y_modes_intervention = y_modes_intervention.median() if is_median else y_modes_intervention.mean()
#             y_modes_post_intervention = y_modes.rolling(n, min_periods=1, center=True)
#             y_modes_post_intervention = y_modes_post_intervention.median() if is_median else y_modes_post_intervention.mean()
#             y_modes_post_intervention = y_modes_post_intervention.loc[cons_post_intervention.index.get_values()]
#             return concat([y_modes_intervention, y_modes_post_intervention])
#
#         for p, n in product(Defs._glucose_adj_neighbors_percentiles, Defs._glucose_adj_neighbors):
#             y_modes = x['GlucoseValue'].groupby('ConnectionID').apply(Defs.find_most_frequent_cgm_value)
#             # y_modes_rolling = y_modes.rolling(n, min_periods=1, center=True).median()
#             y_modes_rolling = rolling(y_modes, n, cons_intervention, cons_post_intervention, is_median=True)
#             gvs.append(x['GlucoseValue'] + y_modes_rolling - y_modes)
#             gv_names.append('GlucoseAdj{}N{}_M'.format(p, n))
#             # y_modes_rolling = y_modes.rolling(n, min_periods=1, center=True).mean()
#             y_modes_rolling = rolling(y_modes, n, cons_intervention, cons_post_intervention, is_median=False)
#             gvs.append(x['GlucoseValue'] + y_modes_rolling - y_modes)
#             gv_names.append('GlucoseAdj{}N{}_Mm'.format(p, n))
#
#             # y_percentiles = x['GlucoseValue'].groupby('ConnectionID').apply(lambda x: np.percentile(x, p))
#             # gvs.append(x['GlucoseValue'] + y_modes.rolling(n, min_periods=1, center=True).median() - y_percentiles)
#             # gv_names.append('GlucoseAdj{}N{}'.format(p, n))
#             # gvs.append(x['GlucoseValue'] + y_modes.rolling(n, min_periods=1, center=True).mean() - y_percentiles)
#             # gv_names.append('GlucoseAdj{}N{}m'.format(p, n))
#
#         for p in Defs._glucose_adj_percentiles:
#             y_percentile = x.groupby('ConnectionID')['GlucoseValue'].apply(lambda x: np.percentile(x, p)).median()
#             gv_p = x.groupby('ConnectionID')['GlucoseValue'].apply(lambda x: self._adjust_connection_cgm_by_percentile(x['GlucoseValue'], y_percentile, p))
#             gvs.append(gv_p)
#             gv_names.append('GlucoseAdj{}'.format(p))
#
#         res = concat(gvs, 1)
#         res.columns = gv_names
#
#         return res
#
#     def _adjust_cgm(self, is_active=1):
#         cgm = self._cgm
#         # if is_active else self._cgm_inactive
#         res = cgm[['GlucoseValue', 'is profiling']].sort_index(level=['UserID', 'GlucoseTimestamp']).groupby('UserID').apply(self._adjust_user_cgm)
#         res = res.join(cgm['GlucoseValue'])
#         res['PPGR'] = res['GlucoseValue'].groupby(level=['UserID', 'ConnectionID']).apply(Defs._ppgr)
#         res['PPGRMin'] = res[Defs.GLUCOSE_ADJ_COL].groupby(level=['UserID', 'ConnectionID']).apply(Defs._ppgr_min)
#         res['PPGRMax'] = res[Defs.GLUCOSE_ADJ_COL].groupby(level=['UserID', 'ConnectionID']).apply(Defs._ppgr_max)
#         return res.join((Series(res.index.get_level_values('GlucoseTimestamp'), index=res.index) - cgm['T0']).apply(lambda x: x.days).to_frame('DaysFromT0'))
#
#     def _get_user_connection_times_from_t0(self, user_id):
#         cgm = self._cgm
#         # if self._randomization.loc[user_id]['IsActive']==1 else self._cgm_inactive
#         con_header = cgm.loc[[user_id]].groupby(['UserID', 'ConnectionID']).head(1).reset_index('GlucoseTimestamp', drop=True).drop('T0', axis='columns').join(self._randomization['T0'])
#         return (con_header['Period_start'] - con_header['T0']).apply(lambda x: x.days)
#
#     def _load_cgm_data(self, is_active=1):
#         res = Load(Defs._cgm_file)
#         res = res.set_index(['UserID', 'ConnectionID', 'GlucoseTimestamp']).sort_index()
#      #   res = res.loc[list(self._randomization[self._randomization['IsActive'] == is_active].index.get_level_values('UserID'))]
#         res = res.loc[res['T0'].dropna().index]
#         res = self._cgm_filter(res)
#         assert res.index.duplicated().sum() == 0
#         return res
#
#     def _cgm_filter(self, cgm):
#         def connection_calendar_days(x):
#             df = x.reset_index('GlucoseTimestamp')
#             return (df['GlucoseTimestamp'][-1] - df['GlucoseTimestamp'][0]).days
#
#         def connection_data_days(x):
#             return (x.shape[0] / 96.)
#
#         print 'Original Connections={}'.format(len(cgm.index.get_level_values('ConnectionID').unique()))
#
#         # cal_days_per_con = cgm['GlucoseValue'].groupby(['UserID', 'ConnectionID']).apply(connection_calendar_days)
#         # filtered_connections = cal_days_per_con[cal_days_per_con >= 5].index.get_level_values('ConnectionID')
#         # cgm = cgm.loc[(slice(None), filtered_connections), :]
#         # print 'PostDays Connections={}'.format(len(cgm.index.get_level_values('ConnectionID').unique()))
#
#         data_days_per_con = cgm['GlucoseValue'].groupby(['UserID', 'ConnectionID']).apply(connection_data_days)
#         filtered_connections = data_days_per_con[data_days_per_con >= 3].index.get_level_values('ConnectionID')
#         cgm = cgm.loc[(slice(None), filtered_connections), :]
#         print 'PostData Connections={}'.format(len(cgm.index.get_level_values('ConnectionID').unique()))
#
#         num_org = cgm.shape[0]
#         cgm = cgm.groupby(['UserID', 'ConnectionID']).apply(lambda x: x.iloc[96:].reset_index(level=['UserID', 'ConnectionID'], drop=True))
#         print 'Removed first day of connections from {} to {}'.format(num_org, cgm.shape[0])
#
#         num_org = cgm.shape[0]
#         cgm = cgm[cgm['GlucoseValue'] > 40]
#         print 'Removed {} min cgm values, from {} to {}'.format(num_org-cgm.shape[0], num_org, cgm.shape[0])
#
#         wrong_prof_t2d_conns = [4116,4448, 5102, 4042, 4229, 4230, 4450, 1607, 1605, 1606, 3001, 3002, 4453, 3509, 3302, 3303, 3508, 3510, 2997, 2137,2229, 3241, 2135,1608,1609, 4541,3511]
#         cgm = cgm[~cgm.index.get_level_values(1).isin(wrong_prof_t2d_conns)]
#         print 'Removed {} wrong profiling values from t2d cross over'.format(len(wrong_prof_t2d_conns))
#         return cgm
#
#     def _external_file_dates(self):
#         df = get_file_dates(Defs._external_files)
#         df.index = [time.strftime('%m/%d/%Y %H:%M:%S', t) for t in df.index.get_level_values('Date')]
#         df.loc[:,'File'] = df['File'].apply(lambda x: os.path.split(x)[1])
#         print df
#
#     def _get_united_randomization_file(self):
#         randomization_prediabetic = pd.read_excel(Defs._randomization_file, sheet_name='Sheet1').rename({'Round': 'BatchNum'}, axis='columns')
#         randomization_prediabetic = randomization_prediabetic[~np.isnan(randomization_prediabetic['UserID'])]
#         randomization_prediabetic['Study'] = 'PreD'
#         randomization_t2d = pd.read_excel(Defs._randomization_file, sheet_name='Sheet2').rename({'Batch': 'BatchNum'}, axis='columns')
#         randomization_t2d = randomization_t2d[~np.isnan(randomization_t2d['UserID'])]
#         randomization_t2d['Study'] = 'T2D'
#         randomization_gdm = pd.read_excel(Defs._randomization_file_gdm, sheet_name = 'full study').rename({'Round': 'BatchNum'}, axis='columns')
#         randomization_gdm = randomization_gdm[~np.isnan(randomization_gdm['UserID'])]
#         randomization_gdm['Study'] = 'GDM'
#         res = concat([randomization_prediabetic, randomization_t2d, randomization_gdm], sort=True).set_index('UserID')
#         res['T0-ShortTerm'] = res[['T0-PersonalizedDiet', 'T0-ControlDiet']].min(axis='columns')
#         res.index = res.index.astype(int)
#         assert res.index.duplicated().sum() == 0
#         return res
#
#     def _load_blood_tests(self):
#         def _handle_duplicate_tests(df, meeting_type):
#             df.loc[df['MeetingTypeID']==meeting_type+0.1, 'MeetingTypeID'] = meeting_type
#             org_df = df[df['MeetingTypeID']!=meeting_type]
#             dup_df = df[df['MeetingTypeID']==meeting_type].sort_values('Answer').groupby(df.index.names).first()
#             return concat([org_df, dup_df]).sort_index()
#
#         res = Defs._read_excel(Defs._blood_tests_file)
#
#         res = res[~res['Answer'].isin((' ', '.'))] # Hack
#         res.loc[:, 'Answer'] = res['Answer'].apply(float)
#         res = res.rename({'userID': 'UserID'}, axis=1).set_index(['UserID', 'Question']).sort_index()
#         res = _handle_duplicate_tests(res, 14)
#         res = _handle_duplicate_tests(res, 5)
#         def convert_dt(x):
#             print x
#             return datetime.datetime.strptime(x.replace('8017', '2017'), '%d.%m.%Y')
#         res.loc[:, 'Timestamp'] = res['Timestamp'].apply(lambda x: convert_dt(x))
#         res = res.reset_index()
#         return res[~res.duplicated()]
#
#     def _get_blood_test(self, blood_test, from_t0, as_series=False):
#         res = self._blood_tests[self._blood_tests['Question'] == blood_test].set_index('UserID')
#         if from_t0:
#             res = res.join(self._user_metadata['T0'])
#             res.loc[:, 'Timestamp'] = (res['Timestamp'] - res['T0']).apply(lambda x: x.days)
#         if as_series:
#             res = res.set_index('Timestamp', append=True)['Answer'].rename(blood_test)
#         return res
#
#     def _cgm_rolling(self, days, averages, adj_neighbors=None):
#         suffix = '' if adj_neighbors is None else '_N{}{}'.format(adj_neighbors, Defs.GLUCOSE_NEIGHBORS_SUFFIX)
#         glucose_col = 'GlucoseValue' if adj_neighbors is None else 'GlucoseAdj50N{}{}'.format(adj_neighbors, Defs.GLUCOSE_NEIGHBORS_SUFFIX)
#         if averages:
#             return self._adjusted_glucose.reset_index('ConnectionID').groupby(level='UserID')[glucose_col]. \
#                 apply(lambda x: x.reset_index('UserID', drop=True).sort_index(). \
#                       resample('15min').mean().  # This fills up the timeline, needed to fill up values per day, otherwise there will be gaps
#                       rolling('{}d'.format(days), min_periods=int(days/3.)*96).mean().\
#                       resample('1d', label='right', closed='right').mean()).rename('{}D{}_avg'.format(days, suffix))
#         else:
#             return self._adjusted_glucose.reset_index('ConnectionID').groupby(level='UserID')[glucose_col]. \
#                 apply(lambda x: x.reset_index('UserID', drop=True).sort_index().\
#                       resample('15min').mean().  # This fills up the timeline, needed to fill up values per day, otherwise there will be gaps
#                       rolling('{}d'.format(days), min_periods=int(days/3.)*96).median().\
#                       resample('1d', label='right', closed='right').median()).rename('{}D{}_med'.format(days, suffix))
#
#     def _prepare_cgm_rolling_averages(self):
#         days = [30, 45, 60, 90]
#         averages = [True]
#         neighbors = [None] + Defs._glucose_adj_neighbors
#         df = concat([self._cgm_rolling(d, a, n) for d,a,n in product(days, averages, neighbors)], 1)
#         df = df.stack().reset_index(['GlucoseTimestamp', None]).join(self._user_metadata['T0']).\
#             rename({'level_2': 'AvgType', 0: 'GlucoseValue'}, axis='columns')
#         df['DaysFromT0'] = (df['GlucoseTimestamp'] - df['T0']).apply(lambda x: x.days)
#         return Defs._change_index(df, self._user_metadata, 'RegistrationCode')
#
#     def _prepare_cgm_rolling_ppgr_percentiles(self):
#         def _user_cgm_rolling(df, days, percentile, glucose_col):
#             if len(df[df['DaysFromT0'] < 0].dropna()) > 0:
#                 cgm_threshold = np.nanpercentile(df[df['DaysFromT0']<0][glucose_col], percentile)
#                 f = lambda x: int(x > cgm_threshold)
#             else:
#                 f = lambda x: np.NaN
#             return df[glucose_col].reset_index('UserID', drop=True).sort_index(). \
#                 apply(f). \
#                 resample('15min').mean().\
#                 rolling('{}d'.format(days), min_periods=int(days / 6.) * 96).mean(). \
#                 resample('1d', label='right', closed='right').mean().rename('{}D_{}p'.format(days, percentile))
#
#         def _cgm_rolling(days, percentile, glucose_col):
#             return self._adjusted_glucose[[glucose_col, 'DaysFromT0']].\
#                 reset_index('ConnectionID', drop=True).groupby(level='UserID').\
#                 apply(lambda x: _user_cgm_rolling(x, days, percentile, glucose_col))
#
#         days = [30]
#         percentiles = [90, 95]
#         glucose_cols = ['PPGR']
#         df = concat([_cgm_rolling(d, t, g) for d,t,g in product(days, percentiles, glucose_cols)], 1)
#         df = df.stack().reset_index(['GlucoseTimestamp', None]).join(self._user_metadata['T0']).\
#             rename({'level_2': 'Percentile', 0: 'GlucoseValue'}, axis='columns')
#         df['DaysFromT0'] = (df['GlucoseTimestamp'] - df['T0']).apply(lambda x: x.days)
#         return Defs._change_index(df, self._user_metadata, 'RegistrationCode')
#
#     def _prepare_cgm_rolling_thresholds(self):
#         def _cgm_rolling(days, threshold, glucose_col):
#             suffix = glucose_col[len('GlucoseValue'):]
#             return self._adjusted_glucose[glucose_col].apply(lambda x: int(x>threshold)).\
#                 reset_index('ConnectionID', drop=True).groupby(level='UserID').\
#                 apply(lambda x: x.reset_index('UserID', drop=True).sort_index().\
#                       resample('15min').mean(). # This fills up the timeline, needed to fill up values per day, otherwise there will be gaps
#                       rolling('{}d'.format(days), min_periods=int(days/6.)*96).mean().\
#                       resample('1d', label='right', closed='right').mean()).rename('{}D{}_{}'.format(days, suffix, threshold))
#
#         days = [10, 30]
#         threshs = [120, 130, 140, 150, 160, 180]
#         glucose_cols = [Defs.GLUCOSE_ADJ_COL]
#         df = concat([_cgm_rolling(d, t, g) for d,t,g in product(days, threshs, glucose_cols)], 1)
#         df = df.stack().reset_index(['GlucoseTimestamp', None]).join(self._user_metadata['T0']).\
#             rename({'level_2': 'Thresh', 0: 'GlucoseValue'}, axis='columns')
#         df['DaysFromT0'] = (df['GlucoseTimestamp'] - df['T0']).apply(lambda x: x.days)
#         return Defs._change_index(df, self._user_metadata, 'RegistrationCode')
#
#     def _prepare_cgm_night_rolling_averages(self):
#         def _cgm_rolling(days, times, glucose_col):
#             gl_in_time = self._adjusted_glucose.reset_index()
#             gl_in_time['Time'] = gl_in_time['GlucoseTimestamp'].apply(lambda x: x.time())
#             gl_in_time = gl_in_time[gl_in_time['Time']<=datetime.time(min(23, times), 0, 0)][gl_in_time['Time']>datetime.time(times-2, 0, 0)]
#             gl_in_time = gl_in_time.set_index(['UserID', 'ConnectionID', 'GlucoseTimestamp'])
#             return gl_in_time.reset_index('ConnectionID').groupby(level='UserID')[glucose_col].\
#                 apply(lambda x: x.reset_index('UserID', drop=True).sort_index().rolling('{}d'.format(days), min_periods = 150).mean().resample('1d', label='right', closed='right').mean()).\
#                 rename('30DN{}_{}_avg'.format(Defs.GLUCOSE_NEIGHBORS, times))
#
#         days = [30]
#         times = [2, 4, 6, 8, 10, 12, 14, 16, 18, 19, 20, 22, 24]
#         glucose_cols = [Defs.GLUCOSE_ADJ_COL]
#         df = concat([_cgm_rolling(d, t, n) for d,t,n in product(days, times, glucose_cols)], 1)
#         df = df.stack().reset_index(['GlucoseTimestamp', None]).join(self._user_metadata['T0']).\
#             rename({'level_2': 'Time', 0: 'GlucoseValue'}, axis='columns')
#         df['DaysFromT0'] = (df['GlucoseTimestamp'] - df['T0']).apply(lambda x: x.days)
#         return Defs._change_index(df, self._user_metadata, 'RegistrationCode')
#
#     def _prepare_user_features(self):
#         features_df1 = self._glucose_features[Defs.GLUCOSE_ADJ_COL]['ByMonth']['ByDay'][['120', '130', '140', '150']].\
#             stack().stack().to_frame('Time above glucose value')
#         features_df2 = self._glucose_features[Defs.GLUCOSE_ADJ_COL]['ByMonth']['ByDay'][['120', '130', '140', '150']].\
#             rename({'120': '120_5N{}'.format(Defs.GLUCOSE_NEIGHBORS_SUFFIX), '130': '130_5N{}'.format(Defs.GLUCOSE_NEIGHBORS_SUFFIX),
#                     '140': '140_5N{}'.format(Defs.GLUCOSE_NEIGHBORS_SUFFIX), '150': '150_5N{}'.format(Defs.GLUCOSE_NEIGHBORS_SUFFIX)}, axis='columns').\
#             stack().stack().to_frame('Time above glucose value')
#         features_df = Defs._change_index(concat([features_df1, features_df2]), self._user_metadata, 'RegistrationCode').reset_index()
#         features_df.rename({'Value': 'Month', 'Threshold': 'GlucoseVal'}, axis='columns', inplace=True)
#         features_df = features_df[features_df['Month']<=30]
#         return features_df
#
#     def _prepare_user_weights(self):
#         weights_df = pd.read_csv(Defs._weight_measures_file).sort_values('MeetingTypeID').set_index('UserID')
#         Defs._meeting_type_2_str(weights_df)
#         weights_df = weights_df[weights_df['MeetingTypeID'] != 'Screening']
#         weights_df = weights_df[weights_df['Weight']>0]
#
#         # Below we fix the first date of T2D to be profiling
#         weights_df = weights_df.join(self._randomization['T0-ShortTerm'])
#         weights_df['DaysFromT0-ShortTerm'] = weights_df['MeetingDate'].astype(np.datetime64) - weights_df['T0-ShortTerm']
#         weights_df.loc[weights_df['DaysFromT0-ShortTerm'].apply(lambda x: x.days) <= 0, 'MeetingTypeID'] = 'Profiling'
#         weights_df = weights_df[weights_df['MeetingTypeID'] != 'PNP1 Connection']
#         weights_df = weights_df.join(self._randomization['Height'])
#         weights_df['BMI'] = weights_df ['Weight'] / weights_df ['Height'] / weights_df ['Height'] * 10000.
#         weights_df = weights_df.reset_index()
#
#         # Hack until we fix the table
#         weights_df.loc[np.logical_or(weights_df['Hips'] < 20, weights_df['Hips'] > 160), 'Hips'] = np.NaN
#         weights_df.loc[np.logical_or(weights_df['Waist'] < 20, weights_df['Waist'] > 160), 'Waist'] = np.NaN
#         weights_df.loc[np.logical_or(weights_df['HeartRate'] < 30, weights_df['HeartRate'] > 120), 'HeartRate'] = np.NaN
#         weights_df.loc[np.logical_or(weights_df['BMR'] == 0, weights_df['BMR'] > 4000), 'BMR'] = np.NaN
#
#         weights_df_table = weights_df[['RegistrationCode', 'Weight', 'TrunkFat_Perc', 'Fat_perc', 'MeetingTypeID', 'BP_sys', 'BP_dia', 'BMR', 'HeartRate', 'Waist', 'Hips', 'BMI']].\
#             groupby(['RegistrationCode', 'MeetingTypeID']).mean().unstack('MeetingTypeID')
#         for f_weight in ['Weight', 'Fat_perc', 'BP_sys', 'BP_dia', 'BMR', 'HeartRate', 'Waist', 'Hips', 'BMI']:
#             weights_df_table[(f_weight, 'Max_Profiling_T0')] = weights_df_table[[(f_weight, 'T0'), (f_weight, 'Profiling')]].max(1)
#             weights_df_table[(f_weight, 'Avg_Profiling_T0')] = weights_df_table[[(f_weight, 'T0'), (f_weight, 'Profiling')]].mean(1)
#         weights_df_table.columns.names = ['WeightFeature', 'MeetingTypeID']
#
#         return weights_df, weights_df_table
#
#     def _fix_blood_tests(self, df):
#         if df.loc[869485, ('HbA1C%', 'T0')] > 9:
#             df.loc[869485, ('HbA1C%', 'T0')] = None
#         if df.loc[479878, ('HbA1C%', 'T0')] > 9:
#             df.loc[479878, ('HbA1C%', 'T0')] = None
#         for blood_test in ['HbA1C%', 'FPG']:
#             df[(blood_test, 'Max_Screening_T0')] = df[[(blood_test, 'T0'), (blood_test, 'Screening')]].max(1)
#             df[(blood_test, 'Avg_Screening_T0')] = df[[(blood_test, 'T0'), (blood_test, 'Screening')]].mean(1)
#
#     def _compute_ppgr_2_time_above(self):
#         thr = 140
#         thr_c = 'TimeAbove{}'.format(thr)
#         adjusted_glucose = self._adjusted_glucose.join(self._user_metadata['ConsPerUser'])
#         adjusted_glucose = adjusted_glucose[adjusted_glucose['ConsPerUser'] >= 6]
#         time_above = adjusted_glucose['GlucoseAdj50'].groupby(level=['UserID', 'ConnectionID']).apply(lambda x: Defs._time_above(x, thr))
#         time_above.name = thr_c
#         time_above = adjusted_glucose[['PPGR', 'GlucoseAdj50']].join(time_above).dropna()
#
#         time_above = time_above[np.logical_and(time_above['GlucoseAdj50'] >= 90, time_above['GlucoseAdj50'] <= 105)]
#         time_above['PPGR'] = time_above['PPGR'].apply(lambda x: 5 * int(x / 5))
#         time_above = time_above.groupby(['PPGR', thr_c]).size().groupby('PPGR').apply(lambda x: x / x.sum()).reset_index()
#         time_above = time_above[time_above['PPGR'] < 100]
#         time_above = time_above.set_index(['PPGR', thr_c]).unstack(thr_c).fillna(0)
#         time_above.columns = time_above.columns.droplevel(0)
#         return time_above
#
#     def _get_regid_order(self):
#         return [int(reg_id) for reg_id in self._user_metadata.sort_values(['T0', 'is main', 'RegistrationCode'])['RegistrationCode'].dropna()]
#
#     def _gen_meals_df(self):
#         user_meals_foods_df = Load(Defs._user_meals_foods_df) \
#             [['Weight', 'Protein_g', 'TotalLipid_g', 'Carbohydrate_g', 'Energy_kcal',
#               'SugarsTotal_g', 'TotalDietaryFiber_g', 'Sodium_mg', 'TotalTransFattyAcids_g', 'TotalSaturatedFattyAcids_g',
#               'TotalMonounsaturatedFattyAcids_g', 'TotalPolyunsaturatedFattyAcids_g',]].\
#             reset_index().\
#             rename({'MealEventID': 'EventID', 'Protein_g': 'Protein', 'TotalLipid_g': 'Lipid', 'Carbohydrate_g': 'Carbs',
#                     'Energy_kcal': 'Energy', 'SugarsTotal_g': 'SugarsTotal',
#                     'TotalDietaryFiber_g': 'Fiber', 'Sodium_mg': 'Sodium', 'TotalTransFattyAcids_g': 'TransFat',
#                     'TotalSaturatedFattyAcids_g': 'SaturatedFat', 'TotalMonounsaturatedFattyAcids_g': 'MonoUnsaturatedFat',
#                     'TotalPolyunsaturatedFattyAcids_g': 'PolyUnsaturatedFat'},
#                    axis='columns').\
#             drop('RunningFoodIndex', axis='columns').\
#             groupby(['UserID', 'EventID']).sum()
#
#         meals_df = Load(Defs._meals_file).reset_index().sort_values('EventType').rename({'RegCode': 'RegistrationCode'}, axis=1)
#         meals_df.rename({'Score max': 'Meal score', 'EventType': 'Meal type'}, axis=1, inplace=True)
#
#         # HACK
#         meals_df['Carbs'] = meals_df['Carbs'].apply(lambda x: eval(x)[0] if type(x) == unicode else x)
#         Defs._meal_type_2_str(meals_df)
#
#         # HACK: Nastya please check why 'Meal score' is not filled out ('Score max' in your original file)
#         meals_df.loc[np.isnan(meals_df['Meal score']), 'Meal score'] = meals_df.loc[np.isnan(meals_df['Meal score']), 'Score_new']
#         meals_df = meals_df.set_index(['UserID', 'EventID'])
#
#         meals_df = meals_df.join(user_meals_foods_df, lsuffix='_')
#
#         # QA: Look for discrepancies in energy reports between the two meal sources
#         energy_diffs = meals_df[((meals_df['Energy_'] - meals_df['Energy']).abs() > 10)]
#         concat([energy_diffs, (energy_diffs['Energy_'] - energy_diffs['Energy']).rename('EnergyDiff')], 1).sort_values('EnergyDiff').\
#             to_excel(os.path.join(Defs._paper_qa_dir, 'energy_large_diffs.xlsx'))
#
#         meals_df = meals_df.drop(['FoodID', 'unit', 'unittype', 'Energy_', 'Protein_', 'Carbs_', 'Lipid_', 'Score', 'Score_new', 'Control'], axis='columns')
#
#         energy_threshold = 2000
#         for nutrient in ['Carbs', 'Protein', 'Lipid', 'Weight', 'SugarsTotal', 'Fiber', 'Sodium', 'TransFat', 'SaturatedFat', 'MonoUnsaturatedFat', 'PolyUnsaturatedFat']:
#             meals_df.loc[meals_df['Energy'] > energy_threshold, nutrient] *= energy_threshold / meals_df.loc[meals_df['Energy'] > energy_threshold, 'Energy']
#         meals_df.loc[meals_df['Energy'] > energy_threshold, 'Energy'] = energy_threshold
#         meals_df = meals_df[~np.isnan(meals_df['Energy'])]
#
#         assert meals_df.index.duplicated().sum()==0
#
#         userid2meal_ppgrs = self._meal_ppgrs.set_index(['UserID', 'EventID'])['Measured_PPRG']
#         assert userid2meal_ppgrs.index.duplicated().sum()==0
#         meals_df = meals_df.join(userid2meal_ppgrs)
#
#         return meals_df
#
#     def _gen_dietary_features(self):
#         res = Defs.DataStruct()
#
#         meal_scores_by_time = self._meals_df.join(self._user_metadata[['T0', 'ProfilingStartByDayFromT0']])
#         meal_scores_by_time['MonthsFromT0'] = (meal_scores_by_time['Timestamp'] - meal_scores_by_time['T0']).apply(lambda x: x.days).apply(Defs._month_groups)
#         meal_scores_by_time['3MonthsFromT0'] = (meal_scores_by_time['Timestamp'] - meal_scores_by_time['T0']).apply(lambda x: x.days).apply(Defs._3month_groups)
#         meal_scores_by_time['6MonthsFromT0'] = (meal_scores_by_time['Timestamp'] - meal_scores_by_time['T0']).apply(lambda x: x.days).apply(Defs._6month_groups)
#         meal_scores_by_time['DaysFromT0'] = (meal_scores_by_time['Timestamp'] - meal_scores_by_time['T0']).apply(lambda x: x.days)
#         meal_scores_by_time['RecentDayGroup'] = (pd.Timestamp(datetime.date.today())-meal_scores_by_time['Timestamp']).apply(lambda x: x.days).apply(lambda x: int(x/float(Defs._num_days_for_recent_day_group))+1)
#         meal_scores_by_time = meal_scores_by_time[np.logical_or(meal_scores_by_time['DaysFromT0'] > -30, meal_scores_by_time['DaysFromT0'] >= meal_scores_by_time['ProfilingStartByDayFromT0'])]
#         num_profiling_logging_days = meal_scores_by_time.set_index('RegistrationCode', append=True).groupby('RegistrationCode')['DaysFromT0'].apply(lambda x: len(np.unique(x[x < 0]))).rename('ProfilingDays')
#
#         def _calc_energy_frac_per_meal_score(meal_scores_by_time, groupby_col, max_val=None):
#             res = meal_scores_by_time.groupby(['RegistrationCode', 'Meal score', groupby_col])['Energy'].sum(). \
#                 groupby(level=['RegistrationCode', groupby_col]).apply(lambda x: x / x.sum())
#             return res[(res.reset_index(groupby_col)[groupby_col]<=max_val).values] if max_val is not None else res
#
#         res.energy_frac_per_meal_score_by_1months = _calc_energy_frac_per_meal_score(meal_scores_by_time, 'MonthsFromT0', 15)
#         res.energy_frac_per_meal_score_by_3months = _calc_energy_frac_per_meal_score(meal_scores_by_time, '3MonthsFromT0', 4)
#         res.energy_frac_per_meal_score_by_6months = _calc_energy_frac_per_meal_score(meal_scores_by_time, '6MonthsFromT0', 1)
#         res.energy_frac_per_meal_score_by_day = _calc_energy_frac_per_meal_score(meal_scores_by_time, 'DaysFromT0')
#         res.energy_frac_per_meal_score_by_recent = _calc_energy_frac_per_meal_score(meal_scores_by_time, 'RecentDayGroup', 2)
#
#         def _calc_energy_per_day_old(meal_scores_by_time, max_days_per_group, groupby_col, max_val=None):
#             res = meal_scores_by_time.groupby(['RegistrationCode', groupby_col])['Energy'].sum()
#             num_days = meal_scores_by_time.groupby(['RegistrationCode', groupby_col]).apply(
#                 lambda x: (datetime.date.today() - x.Timestamp).apply(lambda x: min(max_days_per_group, x.days)).max()).to_frame('NumDays').join(num_profiling_logging_days)
#             if max_days_per_group >= 30:
#                 num_days.loc[(slice(None), 0), 'NumDays'] = num_days.loc[(slice(None), 0), :].min(axis='columns')
#             res = (res / num_days['NumDays']).to_frame('Kcals per day').\
#                 join(self._user_metadata[['RegistrationCode', 'Calories']].dropna().set_index('RegistrationCode')['Calories'].rename('Kcals per day % of target'))
#             # res = (meal_scores_by_time.groupby(['RegistrationCode', groupby_col])['Energy'].sum() / \
#             #        meal_scores_by_time.groupby(['RegistrationCode', groupby_col]). \
#             #        apply(lambda x: (datetime.date.today() - x.Timestamp).apply(lambda x: min(max_days_per_group, x.days)).max())).to_frame('Kcals per day').\
#             #     join(self._user_metadata[['RegistrationCode', 'Calories']].dropna().set_index('RegistrationCode')['Calories'].rename('Kcals per day % of target')). \
#             #     join(num_profiling_logging_days)
#             # if max_days_per_group >= 30:
#             #     res.loc[(slice(None), 0), 'Kcals per day'] = res.loc[(slice(None), 0), 'Kcals per day'] * max_days_per_group / res.loc[(slice(None), 0), 'ProfilingDays']
#             # res = res.drop('ProfilingDays', axis='columns')
#             res.loc[:, 'Kcals per day % of target'] = res['Kcals per day']/res['Kcals per day % of target']
#             return res[(res.reset_index(groupby_col)[groupby_col]<=max_val).values] if max_val is not None else res
#
#         def _calc_energy_per_day(meal_scores_by_time, max_days_per_group, groupby_col, max_val=None, div_by_max=False):
#             if max_val is not None:
#                 meal_scores_by_time = meal_scores_by_time[meal_scores_by_time[groupby_col]<=max_val]
#             res = meal_scores_by_time.groupby(['RegistrationCode', groupby_col])['Energy'].sum()
#             num_days = pd.Series(meal_scores_by_time.groupby(['RegistrationCode', groupby_col]).apply(
#                 lambda x: (pd.Timestamp(datetime.date.today()) - x.Timestamp).apply(lambda x: min(max_days_per_group, x.days)).max())).to_frame('NumDays').join(num_profiling_logging_days)
#             if max_days_per_group >= 30:
#                 num_days.loc[(slice(None), 0), 'NumDays'] = num_days.loc[(slice(None), 0), :].min(axis='columns')
#             res = res / float(max_days_per_group) if div_by_max else res / num_days['NumDays']
#             res = res.to_frame('Kcals per day').join(self._user_metadata[['RegistrationCode', 'Calories']].dropna().set_index('RegistrationCode')['Calories'].rename('Kcals per day % of target'))
#             res.loc[:, 'Kcals per day % of target'] = res['Kcals per day']/res['Kcals per day % of target']
#             return res
#
#         energy_per_day_by_1months = _calc_energy_per_day(meal_scores_by_time, 30, 'MonthsFromT0', 15)
#         energy_per_day_by_3months = _calc_energy_per_day(meal_scores_by_time, 90, '3MonthsFromT0', 4)
#         energy_per_day_by_6months = _calc_energy_per_day(meal_scores_by_time, 180, '6MonthsFromT0', 1)
#         energy_per_day_by_day = _calc_energy_per_day(meal_scores_by_time, 1, 'DaysFromT0')
#         energy_per_day_by_recent = _calc_energy_per_day(meal_scores_by_time, Defs._num_days_for_recent_day_group, 'RecentDayGroup', 2, div_by_max=True)
#
#         nutrients_multipliers = {
#             'Protein': 4,
#             'Lipid': 9,
#             'Carbs': 4,
#             'Energy': 1,
#             'SugarsTotal': 1000,
#             'Fiber': 1000,
#             'Sodium': 1,
#             'TransFat': 9,
#             'SaturatedFat': 9,
#             'MonoUnsaturatedFat': 9,
#             'PolyUnsaturatedFat': 9,
#         }
#
#         nutrients_renames = {
#             'Protein': 'Protein % kcal',
#             'Lipid': 'Lipid % kcal',
#             'Carbs': 'Carbs % kcal',
#             'Energy': 'Energy',
#             'SugarsTotal': 'SugarsTotal g/1000 kcal',
#             'Fiber': 'Fiber g/1000 kcal',
#             'Sodium': 'Sodium',
#             'TransFat': 'TransFat % kcal',
#             'SaturatedFat': 'SaturatedFat % kcal',
#             'MonoUnsaturatedFat': 'MonoUnsaturatedFat % kcal',
#             'PolyUnsaturatedFat': 'PolyUnsaturatedFat % kcal',
#         }
#
#         def _calc_energy_frac_per_nutrient(meal_scores_by_time, groupby_col, max_val=None):
#             nutrients_df = meal_scores_by_time.groupby(['RegistrationCode', groupby_col])[nutrients_multipliers.keys()].sum()
#             for nutrient, multiplier in nutrients_multipliers.items():
#                 if nutrient is not 'Energy':
#                     nutrients_df.loc[:, nutrient] *= multiplier / nutrients_df.loc[:, 'Energy']
#             nutrients_df = nutrients_df.rename(nutrients_renames, axis='columns')
#             return nutrients_df[(nutrients_df.reset_index(groupby_col)[groupby_col]<=max_val).values] if max_val is not None else nutrients_df
#
#         energy_frac_per_nutrient_by_1months = _calc_energy_frac_per_nutrient(meal_scores_by_time, 'MonthsFromT0', 6)
#         energy_frac_per_nutrient_by_3months = _calc_energy_frac_per_nutrient(meal_scores_by_time, '3MonthsFromT0', 4)
#         energy_frac_per_nutrient_by_6months = _calc_energy_frac_per_nutrient(meal_scores_by_time, '6MonthsFromT0', 1)
#         energy_frac_per_nutrient_by_day = _calc_energy_frac_per_nutrient(meal_scores_by_time, 'DaysFromT0')
#         energy_frac_per_nutrient_by_recent = _calc_energy_frac_per_nutrient(meal_scores_by_time, 'RecentDayGroup', 2)
#
#         def _calc_main_control_scores(meal_scores_df, nutrients_df, groupby_col):
#             main_score = meal_scores_df.reset_index(['RegistrationCode', groupby_col]).join(Defs._main_score_table.set_index('Meal score')).dropna().\
#                 groupby(['RegistrationCode', groupby_col]).apply(lambda x: (x['Energy'] * x['Main score']).sum()).rename('main_score')
#
#             control_score_df = nutrients_df[['Carbs % kcal', 'Lipid % kcal', 'SaturatedFat % kcal', 'Fiber g/1000 kcal']]
#             control_score_carbs = Series(np.round(np.interp(100*control_score_df['Carbs % kcal'], [25,45], [0,100])), index=control_score_df['Carbs % kcal'].index, name='control_score_carbs')
#             control_score_lipid = Series(np.round(np.interp(100*control_score_df['Lipid % kcal'], [35,50], [100,0])), index=control_score_df['Lipid % kcal'].index, name='control_score_fat')
#             control_score_saturated_fat = Series(np.round(np.interp(control_score_df['SaturatedFat % kcal'], [0.1,0.2], [100,0])), index=control_score_df['SaturatedFat % kcal'].index, name='control_score_saturated_fat')
#             control_score_fiber = Series(np.round(np.interp(control_score_df['Fiber g/1000 kcal'], [0,15], [0,100])), index=control_score_df['Fiber g/1000 kcal'].index, name='control_score_fiber')
#             control_score_macronutrients = concat([control_score_carbs, control_score_lipid], 1).mean(axis='columns').rename('control_score_macronutrients')
#             control_score = (concat([control_score_carbs*2., control_score_lipid*2., control_score_saturated_fat, control_score_fiber], 1).sum(axis='columns')/6.).rename('control_score')
#
#             return concat([main_score, control_score_carbs, control_score_lipid, control_score_saturated_fat, control_score_fiber, control_score_macronutrients, control_score], 1)
#
#         main_control_scores_by_1months = _calc_main_control_scores(res.energy_frac_per_meal_score_by_1months, energy_frac_per_nutrient_by_1months, 'MonthsFromT0')
#         main_control_scores_by_3months = _calc_main_control_scores(res.energy_frac_per_meal_score_by_3months, energy_frac_per_nutrient_by_3months, '3MonthsFromT0')
#         main_control_scores_by_6months = _calc_main_control_scores(res.energy_frac_per_meal_score_by_6months, energy_frac_per_nutrient_by_6months, '6MonthsFromT0')
#         main_control_scores_by_day = _calc_main_control_scores(res.energy_frac_per_meal_score_by_day, energy_frac_per_nutrient_by_day, 'DaysFromT0')
#         main_control_scores_by_recent = _calc_main_control_scores(res.energy_frac_per_meal_score_by_recent, energy_frac_per_nutrient_by_recent, 'RecentDayGroup')
#
#         res.by_1months = concat([energy_frac_per_nutrient_by_1months, energy_per_day_by_1months, main_control_scores_by_1months], 1)
#         res.by_3months = concat([energy_frac_per_nutrient_by_3months, energy_per_day_by_3months, main_control_scores_by_3months], 1)
#         res.by_6months = concat([energy_frac_per_nutrient_by_6months, energy_per_day_by_6months, main_control_scores_by_6months], 1)
#         res.by_day = concat([energy_frac_per_nutrient_by_day, energy_per_day_by_day, main_control_scores_by_day], 1)
#         res.by_recent = concat([energy_frac_per_nutrient_by_recent, energy_per_day_by_recent, main_control_scores_by_recent], 1)
#
#         return res
#
#     def _prepare_time_above_effect_on_mean(self):
#         df = self._adjusted_glucose
#         profiling_cgm = df[df['DaysFromT0'] <= 0][Defs.GLUCOSE_ADJ_COL]
#         baseline = profiling_cgm.groupby(level='UserID').mean()
#         baseline_size = profiling_cgm.groupby(level='UserID').count()
#         dfs = []
#         for thresh in range(120, 165, 5):
#             df = (1.0 - profiling_cgm[profiling_cgm <= thresh].groupby(level='UserID').mean() / baseline).to_frame('D0')
#             df['Thresh'] = thresh
#             df = df.set_index('Thresh', append=True)
#             dfs.append(df)
#         df = concat(dfs).unstack('Thresh')
#         df.columns = df.columns.droplevel(0)
#         df = df.join(baseline.rename('Mean')).join(baseline_size.rename('NumGlucoseValues'))
#         return df
#
#     def _prepare_user_plot_dfs(self):
#         res = Defs.DataStruct()
#         res.stat_col = 'PPGR'
#
#         res.time_above_effect_on_meal = self._prepare_time_above_effect_on_mean()
#
#         res.cgm_rolling_ppgr_percentiles = self._prepare_cgm_rolling_ppgr_percentiles()
#         res.cgm_rolling_averages = self._prepare_cgm_rolling_averages()
#         res.cgm_rolling_threshs = self._prepare_cgm_rolling_thresholds()
#         # res._cgm_night_rolling_averages = self._prepare_cgm_night_rolling_averages()
#
#         res.hours_df = self._adjusted_glucose[[res.stat_col, 'DaysFromT0']].join(\
#             Series(self._adjusted_glucose[res.stat_col].reset_index('GlucoseTimestamp')['GlucoseTimestamp'].\
#                    apply(lambda x: x.hour).values, index=self._adjusted_glucose.index, name='Hour')).\
#                 set_index('Hour', append=True)
#         res.hours_df['DaysFromT0'] = res.hours_df['DaysFromT0'].apply(Defs._month_groups_str).dropna()
#         res.hours_df = Defs._change_index(res.hours_df, self._user_metadata, 'RegistrationCode').sort_index()
#
#         res.features_df = self._prepare_user_features()
#
#         res.weights_df, res.weights_df_table = self._prepare_user_weights()
#
#         res.blood_tests_df = self._blood_tests[['UserID', 'Question', 'Answer', 'MeetingTypeID']].set_index('UserID')
#         Defs._meeting_type_2_str(res.blood_tests_df)
#         res.blood_tests_df = Defs._change_index(res.blood_tests_df, self._user_metadata, 'RegistrationCode').reset_index().rename({'Question': 'BloodTest'}, axis=1)
#         res.blood_tests_df['Answer'] = res.blood_tests_df['Answer'].apply(lambda x: x if isnumeric(x) else None)
#         res.blood_tests_df = res.blood_tests_df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
#         res.blood_tests_table = res.blood_tests_df[res.blood_tests_df.MeetingTypeID.apply(lambda x: x not in [5.1, 14.1])]\
#             [res.blood_tests_df.BloodTest.apply(lambda x: x not in ['Sodium - blood', 'Potassium - blood', 'Calcium, total - blood'])]
#         res.blood_tests_table = res.blood_tests_table.set_index(['RegistrationCode', 'BloodTest', 'MeetingTypeID']).iloc[:,0].groupby(
#             ['RegistrationCode', 'BloodTest', 'MeetingTypeID']).mean().unstack('BloodTest').unstack('MeetingTypeID')
#         self._fix_blood_tests(res.blood_tests_table)
#
#         res.hba1c_df = res.blood_tests_table['HbA1C%'].stack('MeetingTypeID').reset_index()
#         res.hba1c_df = res.hba1c_df[res.hba1c_df.MeetingTypeID.isin(['Screening', 'T0', '3M', '6M', '9M', '12M'])]
#         res.hba1c_df['MeetingTypeID'] = res.hba1c_df['MeetingTypeID'].astype('category').cat.reorder_categories(['Screening', 'T0', '3M', '6M', '9M', '12M'])
#         res.hba1c_df = res.hba1c_df.sort_values('MeetingTypeID')
#         res.hba1c_df.columns = ['RegistrationCode', 'MeetingTypeID', 'HbA1C%']
#
#         # res.ppgr_2_time_above = self._compute_ppgr_2_time_above()
#         # res.ppgr_boundaries = DataFrame(Load(Defs._boundaries_file)).T
#         # res.ppgr_boundaries.index.names = ['RegistrationCode']
#
#         return res
#
#     def _extract_all_hba1c_dates(self):
#         self._blood_tests[self._blood_tests['Question'] == 'HbA1C%'].Timestamp.value_counts().sort_index().to_excel(
#             '/Users/eran/hba1c_dates.xlsx')
#
# if __name__ == '__main__':
#     data = AnalysisData()