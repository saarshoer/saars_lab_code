import os
import glob
import pickle
import pandas as pd
import xgboost as xgb
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from sklearn.model_selection import GridSearchCV

from anti_mwas_functions import gen_cov_f, gen_y_f
# from anti_mwas import df_dir as train_df_dir
# from anti_mwas_lifeline import df_dir as test_df_dir

run_type = 'within'
models_dir = 'models'

train_run_dir = os.path.join(os.path.dirname(train_df_dir), run_type)
test_run_dir = os.path.join(os.path.dirname(test_df_dir), run_type)


def create_model(train, test, fname=None):

    x_train, y_train = train.iloc[:, :-1], train.iloc[:, -1]

    model = GridSearchCV(xgb.XGBRegressor(random_state=42, n_jobs=1), refit=True, cv=5,
                         param_grid={
                             'n_estimators': [1000],
                             'learning_rate': [0.01, 0.05],
                             'max_depth': [2, 10],
                             'min_child_weight': [100, 150]
                         })
    model.fit(x_train, y_train)
    if fname is not None:
        model.best_estimator_.save_model(os.path.join(train_run_dir, models_dir, fname))

    if test is not None:
        x_test, y_test = test.iloc[:, :-1], test.iloc[:, -1]
        test_score = model.score(x_test, y_test)
    else:
        test_score = None
    return test_score, model.best_score_, model.best_estimator_.feature_importances_


def score_models(sy, sx):
    
    def get_data(run_dir, df_dir):
        file1 = os.path.join(run_dir, 'raw_data', f'mb_gwas_{sx if run_type == "within" else "Rep_all"}_{sy}.h5')
        file2 = os.path.join(run_dir, 'raw_data_not_sig', f'mb_gwas_{sx if run_type == "within" else "Rep_all"}_{sy}.h5')
        files = [file1] if df_dir == train_df_dir else [file1, file2]
        xdata = []
        for file in files:
            if os.path.exists(file):
                xdata.append(pd.read_hdf(file).to_frame('MAF'))
                xdata[-1] = xdata[-1][xdata[-1].index.get_level_values('Species') == sx]  # makes it much more memory efficient
        xdata = pd.concat(xdata)
        xdata['Position'] = xdata.index.get_level_values('Position').astype(int)
        xdata = xdata.reset_index('Position', drop=True).set_index('Position', append=True)
        xdata = xdata['MAF'].unstack('SampleName').T
        sig_df = pd.read_hdf(os.path.join(train_run_dir, 'mb_gwas_significant_validation_clumping.h5'))
        sig_df = sig_df[sig_df.index == sig_df['clumping']]
        xdata = xdata.loc[:, sig_df[(sig_df.index.get_level_values('Y') == sy) &
                                    (sig_df.index.get_level_values('Species') == sx)].index]
        if xdata.shape[1] == 0:
            return None, None
        cdata = gen_cov_f(df_dir, [sx], run_type == 'within').df.dropna()
        ydata = gen_y_f(df_dir, [sx], run_type == 'within').df[sy].dropna()

        data = xdata.join(cdata, how='inner').join(ydata, how='inner').astype(float).sample(frac=1, random_state=42)  # order matters

        return data, cdata.shape[1]

    train_data, n_covariates = get_data(train_run_dir, train_df_dir)
    if train_data is None:
        return
    # test_data = get_data(test_run_dir, test_df_dir)[0][train_data.columns]

    results = {}
    results['validation_base_score'], results['base_score'], _ = create_model(train_data.iloc[:, -n_covariates-1:], None, f'{sx}_{sy}_base.json')
    results['validation_snps_score'], results['snps_score'], results['feature_importance'] = create_model(train_data, None, f'{sx}_{sy}_snps.json')
    results['feature_importance'] = results['feature_importance'][:-n_covariates]
    results['feature_names'] = train_data.columns[:-n_covariates-1]

    results_file = os.path.join(train_run_dir, models_dir, f'{sx}_{sy}.pkl')
    with open(results_file, 'wb') as f:
        pickle.dump(results, f)


if __name__ == "__main__":

    # queue
    jobs_dir = os.path.join(train_run_dir, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(os.path.join(train_run_dir, models_dir))

    with qp(jobname=f'RK{run_type[0]}models', _tryrerun=True, _mem_def='10G') as q:
        q.startpermanentrun()
        tkttores = {}

        sig_df = pd.read_hdf(os.path.join(train_run_dir, 'mb_gwas_significant_validation_clumping.h5')).assign(feature_importance=None)
        sig_df = sig_df[sig_df.index == sig_df['clumping']]
        # sig_df = sig_df[sig_df.index.get_level_values('Species') == 'Rep_595']

        print('start sending jobs')
        for i, (speciesY, speciesX) in sig_df.reset_index()[['Y', 'Species']].drop_duplicates().iterrows():
            if not os.path.exists(os.path.join(train_run_dir, models_dir, f'{speciesX}_{speciesY}.pkl')):
                tkttores[i] = q.method(score_models, (speciesY, speciesX))
        print('finished sending jobs')

        print('start waiting for jobs')
        for k, v in tkttores.items():
            q.waitforresult(v)
        print('finished waiting for jobs')

        print('start df update')
        models = pd.DataFrame(columns=['Y', 'Species', 'base_score', 'snps_score'])#'validation_base_score', 'validation_snps_score',
        results_files = glob.glob(os.path.join(train_run_dir, models_dir, f'*_*.pkl'))
        for r_file in results_files:
            with open(r_file, 'rb') as f:
                r = pickle.load(f)
            sig_df.loc[r['feature_names'], 'feature_importance'] = r['feature_importance']
            new_row = models.shape[0]
            models.loc[new_row, ['Y', 'Species']] = r['feature_names'][0][0], r['feature_names'][0][1]
            # models.loc[new_row, 'validation_base_score'] = r['validation_base_score']
            # models.loc[new_row, 'validation_snps_score'] = r['validation_snps_score']
            models.loc[new_row, 'base_score'] = r['base_score']
            models.loc[new_row, 'snps_score'] = r['snps_score']
        models = models.set_index(['Y', 'Species']).dropna(how='all', axis=1)
        models[['base_score', 'snps_score']] = models[['base_score', 'snps_score']].astype(float)
        models.to_hdf(os.path.join(train_run_dir, f'{models_dir}.h5'), key='snps', complevel=9)
        sig_df.to_hdf(os.path.join(train_run_dir, f'mb_gwas_significant_validation_clumping_{models_dir}.h5'), key='snps', complevel=9)
        print('finished df update')
