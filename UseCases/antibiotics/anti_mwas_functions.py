import os
import pandas as pd
from LabData.DataLoaders.Loader import LoaderData


def gen_cov_f(df_dir, species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'meta.df'))[['age', 'gender']]
    # pca = pd.read_pickle(os.path.join(pca_dir, f'{species[0]}.df'))
    # df = df.join(pca)
    if not within:
        df2 = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))[species]
        # for second type of permutation
        # df2 = pd.read_pickle(os.path.join(os.path.dirname(df_dir), 'abundance.df'))[species]
        df = df.\
            join(df2.rename(columns={species[0]: 'abundance'}))#.\
            # join(df2.rename(columns={species[0]: 'MAF_abundance'}))
    return LoaderData(df, None)


def gen_y_f(df_dir, species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))
    if within:
        df = df[species]
    else:
        if species[0] in df.columns:
            df = df.drop(species[0], axis=1)
    return LoaderData(df, None)


# def is_y_valid(y, max_on_most_freq_val_in_col, min_on_non_freq_val_for_y, y_binary):
#     count_most = y.value_counts().max()
#     return (count_most <= max_on_most_freq_val_in_col * len(y)) & \
#            (count_most >= (1 - max_on_most_freq_val_in_col) * len(y)) & \
#            (len(y) - count_most >= min_on_non_freq_val_for_y)