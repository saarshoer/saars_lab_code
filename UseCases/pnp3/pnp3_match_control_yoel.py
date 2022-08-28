import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import scipy.stats as stats
import os

from LabData.DataLoaders.BloodTestsLoader import BloodTestsLoader
from LabData.DataLoaders.BodyMeasuresLoader import BodyMeasuresLoader
from LabData.DataLoaders.GutMBLoader import GutMBLoader
from LabData.DataLoaders.MetabolonLoader import MetabolonLoader
from LabData.DataMergers.DataMerger import DataMerger



blood_loader = BloodTestsLoader().get_data(study_ids=['PNP1', 'PNP3'])
body_loader = BodyMeasuresLoader().get_data(study_ids=['PNP1', 'PNP3'])
gut_mb_loader = GutMBLoader().get_data('segata_species',take_log=True, study_ids=['PNP1', 'PNP3'], min_col_present_frac=0.1)
serume_metabolomics_loader = MetabolonLoader().get_data(metabolon_runs=['Metabolon_MAR2017', 'Metabolon_APR2020'], study_ids=['PNP1', 'PNP3'],
robust_zs=True, clip_outliers_std=8, clip_outliers_or_na='clip', fill_missing_with_min=True,
log_transform_base=10, norm_within_run_using_anchors=True, norm_across_run_using_anchors=True)

blood_data = blood_loader.df
body_data = body_loader.df
gut_mb_data = gut_mb_loader.df
metabloites_data = serume_metabolomics_loader.df
blood_meta = blood_loader.df_metadata
body_meta = body_loader.df_metadata
gut_mb_meta = gut_mb_loader.df_metadata
metabloites_meta = serume_metabolomics_loader.df_metadata
cytokines_data = pd.read_csv(
    '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/data_frames/cytokines/cytokines-nastya.csv')
cytokines_data.reset_index(inplace=True)
cytokines_data = cytokines_data.astype({'RegistrationCode':str},errors='raise')
cytokines_data.set_index(['Plate ID'], inplace=True)
cytokines_meta = cytokines_data.copy()
carbs = pd.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/data_frames/diet_20_corrected.df')[['% Carbohydrates']]
carbs.reset_index(inplace=True)
carbs_6 = carbs[carbs['time_point'] == '6months']
carbs_6.set_index(['person'], inplace=True)



#syncs data 1 to data 2
def sync_data(data1, data2):
    data1.reset_index(inplace=True)
    data2.reset_index(inplace=True)
    data1.set_index(['RegistrationCode'], inplace=True)
    data2.set_index(['RegistrationCode'], inplace=True)
    for i in data1.index:
        if not (i in data2.index):
            data1.drop(index=i, inplace=True)
    return data1

#return first test date for a given person
def find_first_date(reg_code, data):
    dates = []
    for i in data.index:
        if reg_code == i[0]:
            dates.append(i[1])
    return min(dates)

#returns the data with base tests only
def first_test_only(data):
    original_index = []
    for i in data.index:
        if not (i[0] in original_index):
            original_index.append(i[0])
    first_test_dates = [find_first_date(person, data) for person in original_index]
    first_index = zip(original_index, first_test_dates)
    return data.loc[first_index].reset_index()

#creates the sub cohort index list
def create_sub_dataframe(PNP1_blood_meta, PNP3_blood_meta, BMI_flag, person, body_data
                         , age_diff=3, bmi_diff=5):
    if len(PNP3_blood_meta.loc[person]) == 2: #sometimes there are two results for the same test
        person_gender = PNP3_blood_meta.loc[person, 'gender'][0]
        person_age = PNP3_blood_meta.loc[person, 'age'][0]
    else:
        person_gender = PNP3_blood_meta.loc[person, 'gender']
        person_age = PNP3_blood_meta.loc[person, 'age']
    # parsing by age and gender
    sub_df = PNP1_blood_meta[(PNP1_blood_meta['gender'] == person_gender) &
                             (np.abs(PNP1_blood_meta['age'] - person_age) <= age_diff)]
    if BMI_flag: #parse by BMI
        person_BMI = body_data.loc[person, 'bmi']
        sub_body_df = body_data.loc[sub_df.index]
        sub_body_df = sub_body_df[np.abs(sub_body_df['bmi'] - person_BMI) <= bmi_diff]
        sub_df = sub_df.loc[sub_body_df.index]
    return sub_df.index.to_series()

#creates a data frame with sub cohort index for each PNP3 person

#each row is a PNP3 person, first colomn is the sub cohort list,
# second is the sub cohort length

#DROPS: PNP3 persons how couldnt match with at least 3 PNP1 persons

def create_sub_df_dictionary(blood_data, blood_meta, body_data, bmi_flag, age_diff=5, bmi_diff=3):
    base_blood_meta = first_test_only(blood_meta).set_index(['RegistrationCode'])
    base_body_data = first_test_only(body_data).set_index(['RegistrationCode'])
    base_blood_data = first_test_only(blood_data).set_index(['RegistrationCode'])
    base_blood_meta = sync_data(base_blood_meta, base_body_data)#sync body and blood
    base_body_data = sync_data(base_body_data, base_blood_meta)#sync body and blood
    base_blood_data.reset_index(inplace=True)
    base_blood_data.set_index(['RegistrationCode'], inplace=True)
    PNP1_blood_data = base_blood_data.loc[base_blood_meta[base_blood_meta['StudyTypeID'] == 1].index]
    PNP1_blood_data = PNP1_blood_data[PNP1_blood_data['bt__hba1c'] < 5.7]#parsing by non diabitic
    PNP1_blood_meta = base_blood_meta.loc[PNP1_blood_data.index]
    # parsing by base PNP3 test - which is research stage 5
    PNP3_blood_meta = blood_meta[(blood_meta['StudyTypeID'] == 3) & (blood_meta['research_stage'] == 5)].reset_index()
    PNP3_blood_meta.set_index(['RegistrationCode'], inplace=True)
    dictionary = [create_sub_dataframe(PNP1_blood_meta, PNP3_blood_meta, bmi_flag, person, base_body_data) for person in
                  PNP3_blood_meta.index] #creating sub cohort for each person
    length = [len(l) for l in dictionary] #for sub cohort length histogram
    dictionary_df = pd.DataFrame({'ReferenceIndex': dictionary, 'ReferenceLength': length}, index=PNP3_blood_meta.index)
    dictionary_df = dictionary_df[dictionary_df['ReferenceLength'] >= 3]
    return dictionary_df


#returns the normalized deviations list for a given person and data(blood or body)
def find_blood_body_deviations(dictionary, data, meta, person, stage):
    print(person)#to feel the time passes by...
    deviations = []
    PNP1 = data.loc[meta[meta['StudyTypeID'] == 1].index]
    PNP1 = first_test_only(PNP1).set_index(['RegistrationCode']).drop(columns='Date')
    #dropping date becasue its a numerical column that we dont want
    PNP1.reset_index(inplace=True)
    PNP1.set_index('RegistrationCode', inplace=True)
    #parsing by current stage and study type
    PNP3 = data.loc[meta[(meta['research_stage'] == stage) & (meta['StudyTypeID'] == 3)].index]
    PNP3.reset_index(inplace=True)
    PNP3.set_index('RegistrationCode', inplace=True)
    if not (person in PNP3.index): #if that person did not test in the current reasearch stage
        return np.nan
    sub_data = PNP1.loc[dictionary.loc[person, 'ReferenceIndex']] #sub cohort data
    std_list = sub_data.std()#std list
    mean_list = sub_data.mean()#mean list
    for feature in std_list.index:
        if type(PNP3.loc[person, feature]) == np.float64 or type(PNP3.loc[person, feature]) == float:
            deviations.append((PNP3.loc[person, feature] - mean_list[feature]) / std_list[feature])#normalizing
        else:
            deviations.append((PNP3.loc[person, feature].mean() - mean_list[feature]) / std_list[feature])#normalizing, if machine returned 2 results for the same test
    return pd.Series(deviations, index=std_list.index)

#returns the normalized deviations list for a given person and data(mb,mata)
def find_mb_metabolon_cytokines_dev(dic, data, meta, person, stage, dtype_str):
    print(person)#to feel the time passes by...
    if stage == 9:#no test for this stage
        return np.nan
    if dtype_str == 1:#determeine the data type
        sample_index = 'SerumName'
    if dtype_str == 0:
        sample_index = 'SampleName'
    if dtype_str == 2:
        sample_index = 'Plate ID'
    deviations = []
    PNP3_meta = meta[meta['StudyTypeID'] == 3].reset_index()
    PNP3_meta = PNP3_meta.set_index(['RegistrationCode', 'StorageDT'])
    PNP3_meta_first = first_test_only(PNP3_meta).set_index(['RegistrationCode', 'StorageDT'])#base test
    PNP3_meta_last = PNP3_meta[~PNP3_meta.index.isin(PNP3_meta_first.index)]#6 mothes test
    if stage == 14:#parse by stage
        PNP3_meta = PNP3_meta_last
    if stage == 5:#parse by stage
        PNP3_meta = PNP3_meta_first
    PNP3_meta.reset_index(inplace=True)
    PNP3_meta.set_index(['RegistrationCode'], inplace=True)
    if not (person in PNP3_meta.index):#meaning he did not test in current stage
        return np.nan
    PNP1_meta = meta[meta['StudyTypeID'] == 1].reset_index().drop_duplicates(subset=['RegistrationCode'])
    PNP1_meta.set_index(['RegistrationCode'], inplace=True)
    cohort_index = [j for j in dic.loc[person, 'ReferenceIndex'] if j in PNP1_meta.index]#refernce cohort, but only those who have done this test
    cohort_sample_index = [PNP1_meta.loc[j, sample_index] for j in cohort_index]
    if len(cohort_sample_index) < 5:#reference cohort too small
        return np.nan
    sub_data = data.loc[cohort_sample_index]
    std_list = sub_data.std()#std list
    mean_list = sub_data.mean()#mean list
    for feature in std_list.index:
        if type(data.loc[PNP3_meta.loc[person, sample_index], feature]) == np.float64 or type(
                data.loc[PNP3_meta.loc[person, sample_index], feature]) == float:
            deviations.append((data.loc[PNP3_meta.loc[person, sample_index], feature] - mean_list[feature]) / std_list[feature])#normalizing
        else:
            deviations.append((data.loc[PNP3_meta.loc[person, sample_index], feature].mean() - mean_list[feature]) / std_list[feature])#normalizing, if machine returned 2 results for the same test
    return pd.Series(deviations, index=std_list.index)

#creats the normalized deviations data frame for a given stage
#each PNP3 person has his daviations list for each meassure(blood, body...)
def stage_deviation_directory(blood_data, blood_meta, body_data, body_meta, gut_mb_data, gut_mb_meta, metabloites_data,
                              metabloites_meta, cytokines_data, cytokines_meta, sub_df_dictionary, research_stage):
    #blood deviations list for each PNP3 person


    blood_deviations = pd.Series(
        [find_blood_body_deviations(sub_df_dictionary, blood_data, blood_meta, i, research_stage)
         for i in sub_df_dictionary.index], index=sub_df_dictionary.index)
    # body deviations list for each PNP3 person
    body_deviations = pd.Series([find_blood_body_deviations(sub_df_dictionary, body_data, body_meta, i, research_stage)
                                 for i in sub_df_dictionary.index], index=sub_df_dictionary.index)
    # mb deviations list for each PNP3 person
    gut_mb_deviations = pd.Series(
        [find_mb_metabolon_cytokines_dev(sub_df_dictionary, gut_mb_data, gut_mb_meta, i, research_stage, 0)
         for i in sub_df_dictionary.index], index=sub_df_dictionary.index)
    # metabloites deviations list for each PNP3 person

    metabloites_deviations = pd.Series(
        [find_mb_metabolon_cytokines_dev(sub_df_dictionary, metabloites_data, metabloites_meta, i, research_stage, 1)
         for i in sub_df_dictionary.index], index=sub_df_dictionary.index)

    # cytokines deviations list for each PNP3 person
    cytokines_deviations = pd.Series(
        [find_mb_metabolon_cytokines_dev(sub_df_dictionary, cytokines_data, cytokines_meta, i, research_stage, 2)
         for i in sub_df_dictionary.index], index=sub_df_dictionary.index)
    return pd.DataFrame({'ReferenceIndex': sub_df_dictionary['ReferenceIndex'], 'BloodDeviations': blood_deviations,
                         'BodyDeviations': body_deviations,
                         'Gut_MB_Deviations': gut_mb_deviations, 'MetabloitesDeviations': metabloites_deviations
                        ,'CytokinesDeviations': cytokines_deviations},index=sub_df_dictionary.index)

#creates a deviation data frame for each pnp3 stage
def create_deviation_directory(blood_data, blood_meta, body_data, body_meta, gut_mb_data, gut_mb_meta,
                               metabloites_data, metabloites_meta,cytokines_data,cytokines_meta, sub_df_dictionary):
    stages = [5, 9, 14]
    stages_deviations = []
    for i in stages:
        stages_deviations.append(
            stage_deviation_directory(blood_data, blood_meta, body_data, body_meta, gut_mb_data, gut_mb_meta,
                                      metabloites_data, metabloites_meta,cytokines_data,cytokines_meta, sub_df_dictionary, i))
    return stages_deviations

#saves the deviation data frame list
def save_deviation(bmi_flag,age_diff,bmi_diff):
    dictionary = create_sub_df_dictionary(blood_data, blood_meta, body_data, bmi_flag=bmi_flag,age_diff=age_diff,bmi_diff=bmi_diff).reset_index().drop_duplicates(
        subset=['RegistrationCode'])
    dictionary.set_index(['RegistrationCode'], inplace=True)
    dev = create_deviation_directory(blood_data, blood_meta, body_data, body_meta, gut_mb_data, gut_mb_meta,metabloites_data, metabloites_meta,cytokines_data,cytokines_meta, dictionary)
    for stage in range(len(dev)):
        dev[stage].to_hdf('Stage' + str(stage) + 'Age' + str(age_diff) + 'Bmi' + str(bmi_diff) + '.h5', key='df', mode='w')




#generate heat maps data
def create_heat_map(stages, data_str, stage):
    data_heat_map = []
    data_index = stages[0].loc[stages[0].index[0], data_str].index
    for i in data_index:
        curr_fe = []
        for p in stages[0].index:
            if type(stages[stage].loc[p, data_str]) == float or type(stages[stage].loc[p, data_str]) == np.float64:
                curr_fe.append(np.nan)
            else:
                curr_fe.append(stages[stage].loc[p, data_str][i])
        data_heat_map.append(curr_fe)
    heat_df = pd.DataFrame(data_heat_map, index=data_index, columns=stages[0].index)
    data_heat_map = np.array(data_heat_map)[~np.isnan(data_heat_map)]
    data_heat_map = data_heat_map.flatten()[data_heat_map != np.inf]
    percentile = np.percentile(data_heat_map,95)
    #plt.hist(data_heat_map,bins=range(-10,10))
    #plt.show()
    for j in heat_df.columns:
        isnan = 1
        for i in heat_df.index:
            if not np.isnan(heat_df.loc[i, j]):
                isnan = 0
        if isnan:
            heat_df.drop(j, inplace=True, axis=1)

    drop = []
    for i in heat_df.index:
        count = 0
        length = len(heat_df.columns)
        for j in heat_df.columns:
            if np.isnan(heat_df.loc[i, j]):
                count += 1

        if count / length > 0.15:
            drop.append(i)
    for i in drop:
        heat_df.drop(i, inplace=True)
    #median = heat_df.median(axis=1)
    #print(median)
    #for i in median:
    return heat_df, percentile

#generate abs val data
def abs_value_data(base, after):
    abs_data = []
    for i in after.index:
        curr_f = []
        for p in after.columns:
            if not i in base.index or not p in base.columns:
                curr_f.append(np.nan)
            else:
                stage_0 = base.loc[i, p]
                stage_6 = after.loc[i, p]
                if not np.isnan(stage_6) and not np.isnan(stage_0):
                    curr_f.append(np.abs(stage_6) - np.abs(stage_0))
                else:
                    curr_f.append(np.nan)
        abs_data.append(curr_f)

    return pd.DataFrame(abs_data, index=after.index, columns=after.columns)
#fills rows with median(for cluster)
def fill_row_median(abs_value_diff_data):
    df = abs_value_diff_data.copy()
    diff_median = abs_value_diff_data.median(axis=1)
    for p in df.columns:
        for f in df.index:
            if np.isnan(df.loc[f, p]):
                df.loc[f, p] = diff_median[f]
    return df

#draws the heat maps as we disscused, takes the deviations data frame list as argument
def draw_heat_maps(stages):
    data_str = ['BloodDeviations', 'BodyDeviations', 'Gut_MB_Deviations', 'MetabloitesDeviations','CytokinesDeviations']
    abs_list = []
    for i in data_str:
        print(i)
        sns.set()
        print('0')
        data_heat_map_base = list(create_heat_map(stages, i, 0))
        print('2')
        data_heat_map_6 = list(create_heat_map(stages, i, 2))
        for f in data_heat_map_6[0].index:
            if f not in data_heat_map_base[0].index:
                data_heat_map_6[0].drop(f, inplace=True)
        for f in data_heat_map_base[0].index:
            if f not in data_heat_map_6[0].index:
                data_heat_map_base[0].drop(f, inplace=True)
        for p in data_heat_map_base[0].columns:
            if p not in data_heat_map_6[0].columns:
                data_heat_map_base[0].drop(p, inplace=True, axis=1)
        abs_value_diff_data = abs_value_data(data_heat_map_base[0], data_heat_map_6[0])
        abs_list.append(abs_value_diff_data)
        carbs_pre = []
        for p in abs_value_diff_data.columns:
            if p in carbs_6.index:
                carbs_pre.append((0,carbs_6.loc[p,'% Carbohydrates'],0,0))
            else:
                carbs_pre.append((1,0,0,0))
        carbs_pre = pd.Series(carbs_pre,index=abs_value_diff_data.columns)
        f, (ax1, ax2) = plt.subplots(1,2,sharey='all')
        f.suptitle(i + ' Heat Maps')
        g3 = sns.clustermap(fill_row_median(abs_value_diff_data), cmap="vlag", vmin=-3,
                         vmax=3,col_colors=carbs_pre,
                            mask=abs_value_diff_data.isnull())
        data_heat_map_base[0] = data_heat_map_base[0].reindex(g3.data2d.index)
        data_heat_map_base[0] = data_heat_map_base[0].reindex(g3.data2d.columns,axis='columns')
        data_heat_map_6[0] = data_heat_map_6[0].reindex(g3.data2d.index)
        data_heat_map_6[0] = data_heat_map_6[0].reindex(g3.data2d.columns, axis='columns')
        carbs_pre = carbs_pre.reindex(g3.data2d.columns)
        g1 = sns.heatmap(data_heat_map_base[0], cmap="vlag", vmin=-data_heat_map_base[1],
                         vmax=data_heat_map_base[1], mask=data_heat_map_base[0].isnull()
                         ,yticklabels=False, ax=ax1)
        g1.set_title('base num of std')
        g1.set_facecolor('black')
        g2 = sns.heatmap(data_heat_map_6[0], cmap="vlag", vmin=-data_heat_map_base[1],
                         vmax=data_heat_map_base[1], mask=data_heat_map_6[0].isnull()
                         ,yticklabels=False, ax=ax2)
        g2.set_title('6 monthes num of std')
        g2.set_facecolor('black')
        g3.fig.suptitle('abs value diff '+i)
        plt.show()
    return abs_list
#unfinished work
def find_carbs_diff_cor(abs_list,carbs):
    carbs_base = carbs[carbs['time_point']=='0months']
    carbs_base.set_index(['person'],inplace=True)
    correlations = []
    for i in abs_list:
        i = i.T
        index = [p for p in i.index if p in carbs_base and p in carbs_6]
        i = i.loc[index]
        i['% Carbohydrates diff'] = [carbs_6.loc[p,'% Carbohydrates']-carbs_base.loc[p,'% Carbohydrates']
                                     for p in i.index]
        corr = i.corr()
        p_values = []
        p_values = p_values[p_values['% Carbohydrates diff'] < 0.05]
        corr = corr.loc[p_values.index]
        corr_df = pd.DataFrame({'c_val':corr['% Carbohydrates diff'],'p_val':p_values['% Carbohydrates diff']},index=corr.index)
        print(corr_df.to_string())
        correlations.append(corr_df)
    return correlations



save_deviation(1,6,4)
