import os
import pandas as pd
from LabData.DataLoaders.Loader import LoaderData


def gen_cov_f(df_dir, species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'meta.df'))[['age', 'gender']].assign(coverage=1)

    pca = pd.read_hdf(os.path.join(os.path.dirname(df_dir), 'pcs_covariate', f'mb_snp_pcs_{species[0]}.h5'))
    pca.columns = [f'PC{int(i)+1}' for i in pca.columns]
    df = df.join(pca)

    if not within:
        df2 = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))[species]
        df = df.join(df2.rename(columns={species[0]: 'abundance'}))
    return LoaderData(df, None)


def gen_y_f(df_dir, species, within):
    df = pd.read_pickle(os.path.join(df_dir, 'abundance.df'))
    if within:
        df = df[species]
    else:
        if species[0] in df.columns:
            df = df.drop(species[0], axis=1)
    return LoaderData(df, None)
