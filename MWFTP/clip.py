"""Provides functionality for clipping extreme values.

Based on NormDistCapping.fit_transform().
"""
import pandas as pd


def clip(data: pd.Series,
         window_frac: float = 0.98,
         stds_clip: float = 6,
         stds_remove: float = 10) -> pd.DataFrame:
    """Returns a copy of the series with the values of column removed or
    clipped, according to their distance from the mean.

    Values stds_remove or more away from the mean are removed. Values stds_clip or
    more away from mean are set to mean+=stds_clip."""

    assert 0 < window_frac <= 1, f'{window_frac}'

    # Find the densest window in the series (lowest span).
    # We consider only the values within this window for calculating mean and std,
    # to avoid extreme values.
    series = data.sort_values()
    series = series[~series.isna()]
    window_size = round(len(series) * window_frac)
    windows = series.rolling(window_size)
    spans = windows.max() - windows.min()
    spans = spans[window_size - 1:]  # Remove prefix of NAs.
    argmin = spans.values.argmin()
    assert argmin >= 0, f'argmin={argmin}'
    series = series[argmin:argmin + window_size]

    # Filter according to std of the dense window.
    mean, std = series.mean(), series.std()
    keep_cond = ((data <= mean + std * stds_remove) &
                 (data >= mean - std * stds_remove))
    data = data[keep_cond][:]

    if len(data) == 0:
        raise ValueError('all the values were clipped')

    # Clip values that are more extreme than stds_clip.
    clip_upper = mean + std * stds_clip
    clip_lower = mean - std * stds_clip
    data.loc[data > clip_upper] = clip_upper
    data.loc[data < clip_lower] = clip_lower

    return data
