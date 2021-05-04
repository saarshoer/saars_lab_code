"""An XY function for regression models."""
import statsmodels.api as sm


def regression(x, y):
    """An XY function for regression models."""
    y_binary = len(set(y)) == 2

    if y_binary:
        sm_result = sm.Logit(y, x).fit(method='nm')  # TODO: understand why nm is the best method
    else:
        sm_result = sm.OLS(y, x).fit()

    result = {}
    coef_interval = sm_result.conf_int()

    result['rsquared'] = sm_result.prsquared if y_binary else sm_result.rsquared
    for col in x.columns:
        result[f'{col}_coef'] =  sm_result.params[col]
        result[f'{col}_pval'] = sm_result.pvalues[col]
        result[f'{col}_coef_025'] = coef_interval.loc[col, 0]
        result[f'{col}_coef_975'] = coef_interval.loc[col, 1]

    return result
