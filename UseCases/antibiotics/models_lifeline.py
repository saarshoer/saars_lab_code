import os
import glob
import pickle
import pandas as pd
import xgboost as xgb
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from sklearn.model_selection import GridSearchCV
from anti_mwas_lifeline import gen_cov_f, gen_y_f, df_dir

run_type = 'between'

run_dir = os.path.join(os.path.dirname(df_dir), run_type)


def create_model(X, y, fname=None):
    model = GridSearchCV(xgb.XGBRegressor(random_state=42, n_jobs=1), refit=True, cv=5,
                         param_grid={})#{'max_depth': [4, 6], 'n_estimators': [50, 100]})
    model.fit(X, y)
    if fname is not None:
        model.best_estimator_.save_model(os.path.join(run_dir, 'models', fname))
    return model.best_score_, model.best_estimator_.feature_importances_


def score_models(sy, sx):
    # data = pd.read_hdf(os.path.join(run_dir, 'raw_data', f'mb_gwas_{sx if run_type == "within" else "Rep_all"}_{sy}.h5'))
    # data = data.reset_index().drop_duplicates().set_index(data.index.names).iloc[:, -1]
    # data.to_hdf(os.path.join(run_dir, 'raw_data', f'mb_gwas_{sx if run_type == "within" else "Rep_all"}_{sy}_new.h5'), key='snps', complevel=9)
    # del data
    # return

    xdata = pd.read_hdf(os.path.join(run_dir, 'raw_data', f'mb_gwas_{sx if run_type == "within" else "Rep_all"}_{sy}.h5')).to_frame('MAF')
    xdata = xdata[xdata.index.get_level_values('Species') == sx]  # makes it much more memory efficient
    xdata['Position'] = xdata.index.get_level_values('Position').astype(int)
    xdata = xdata.reset_index('Position', drop=True).set_index('Position', append=True)
    xdata = xdata['MAF'].unstack('SampleName').T
    cdata = gen_cov_f([sx], run_type == 'within').df.dropna()
    ydata = gen_y_f([sx], run_type == 'within').df[sy].dropna()

    data = xdata.join(cdata, how='inner').join(ydata, how='inner').astype(float).sample(frac=1, random_state=42)  # order matters
    n_covariates = cdata.shape[1]
    del xdata, cdata, ydata

    results = {}
    results['base_score'], _ = create_model(data.iloc[:, -n_covariates-1:-1], data.iloc[:, -1], f'{sx}_{sy}_base.json')
    results['snps_score'], results['feature_importance'] = create_model(data.iloc[:, :-1], data.iloc[:, -1], f'{sx}_{sy}_snps.json')
    results['feature_importance'] = results['feature_importance'][:-n_covariates]
    results['feature_names'] = data.columns[:-n_covariates-1]

    results_file = os.path.join(run_dir, 'models', f'{sx}_{sy}.pkl')
    with open(results_file, 'wb') as f:
        pickle.dump(results, f)

    return results_file


if __name__ == "__main__":

    # queue
    jobs_dir = os.path.join(run_dir, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(os.path.join(run_dir, 'models'))

    with qp(jobname='Lmodels', _tryrerun=True, _mem_def='10G') as q:
        q.startpermanentrun()
        tkttores = {}

        sig_df = pd.read_hdf(os.path.join(run_dir, 'mb_gwas_significant_validation.h5')).assign(feature_importance=None)
        # sig_df = sig_df[sig_df.index.get_level_values('Species') == 'Rep_595']
        models = pd.DataFrame(columns=['Y', 'Species', 'base_score', 'snps_score'])

        # all_probs = []
        print('start sending jobs')
        for i, (speciesY, speciesX) in sig_df.reset_index()[['Y', 'Species']].drop_duplicates().iterrows():
            # if not os.path.exists(os.path.join(run_dir, 'models', f'{speciesX}_{speciesY}.pkl')):
            #     print(speciesX, speciesY)
                # if speciesY not in all_probs:
                #     all_probs.append(speciesY)
            tkttores[i] = q.method(score_models, (speciesY, speciesX))
        print('finished sending jobs')

        print('start waiting for jobs')
        results_files = []
        for k, v in tkttores.items():
            results_files.append(q.waitforresult(v))
        print('finished waiting for jobs')
        results_files = glob.glob(os.path.join(run_dir, 'models', f'*_*.pkl'))

        print('start df update')
        for r_file in results_files:
            with open(r_file, 'rb') as f:
                r = pickle.load(f)
            sig_df.loc[r['feature_names'], 'feature_importance'] = r['feature_importance']
            new_row = models.shape[0]
            models.loc[new_row, ['Y', 'Species']] = r['feature_names'][0][0], r['feature_names'][0][1]
            models.loc[new_row, 'base_score'] = r['base_score']
            models.loc[new_row, 'snps_score'] = r['snps_score']
        models = models.set_index(['Y', 'Species'])
        models[['base_score', 'snps_score']] = models[['base_score', 'snps_score']].astype(float)
        models.to_hdf(os.path.join(run_dir, 'models.h5'), key='snps', complevel=9)
        sig_df.to_hdf(os.path.join(run_dir, 'mb_gwas_significant_validation_models.h5'), key='snps', complevel=9)
        print('finished df update')
