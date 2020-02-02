# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : util
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST15r01
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""util.

Routine Listings
----------------
quadrature
mad_with_extrema
mad_adj_by_std
adaptive_binning

"""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# GENERAL
import numpy as np

from scipy.stats import binned_statistic as binned_stats

from astropy import units as u
from astropy.stats import median_absolute_deviation


#############################################################################
# CODE


def mad_with_extrema(x, mnm=None, mxm=None, return_func=False):
    """Calculate Bounded Median Absolute Deviation.

    with a minimum and maximum allowed mad

    Parameters
    ----------
    x: array
        the data on which to calculate the MAD
    mnm: float
        lower bound
        if MAD(x) < mnm: return mnm
    mxm: float
        upper bound
        if MAD(x) > mxm: return mxm
    return_func: bool
        if True, returns a single-parameter function with set mnm & mxm

    Returns
    -------
    MAD: float
        if return is True
        `mad_with_extrema(., mnm=mnm, mxm=mxm, return_func=False)`

    """
    if return_func is True:
        return lambda x: mad_with_extrema(
            x, mnm=mnm, mxm=mxm, return_func=False
        )

    # astropy MAD
    mad = np.nan_to_num(median_absolute_deviation(x))

    if mnm is not None:

        if issubclass(x.__class__, u.Quantity):
            mnm = mnm * x.unit

        try:  # array
            mad[mad < mnm] = mnm
        except TypeError:  # single value
            if mad < mnm:
                mad = mnm

    if mxm is not None:

        if issubclass(x.__class__, u.Quantity):
            mxm = mxm * x.unit

        try:
            mad[mad > mxm] = mxm
        except TypeError:
            if mad > mxm:
                mad = mxm

    return mad


# /def


def mad_adj_by_std(bin_mad, x, binnumber):
    """Adjust bin_mad >= max(std(x(bin)))."""
    for n in set(binnumber):
        bin_mad[n - 1] = max(np.mean(x[binnumber == n]).value, bin_mad[n - 1])
    return bin_mad


# /def


# ----------------------------------------------------------------------------


def adaptive_binning(
    x,
    y,
    statistic="mean",
    ibins=10,
    pcttol=10,
    minpnts=25,
    xunits=None,
    yunits=None,
):
    """Adaptive binning.

    Parameters
    ----------
    x: array
        the x data
    y: array
        the data to be binned
    statistic: str or func
        statistic method for scipy.stats.binned_statistic
        see binned_statistics documentation for options
    pcttol: float
        percent difference tolerance to stop binning
    minpnts: int
        minimum number of points in a bin

    Returns
    -------
    bin_stat: ndarray
    bin_cents: ndarray
    bin_edges: ndarray
    binnumber: ndarray
    pctdiff: ndarray

    """
    _stats = binned_stats(x, y, statistic="median", bins=ibins)
    bin_stat, bin_edges, binnumber = _stats

    # percent difference between medians
    pctdiff = np.array([0, *np.abs(np.diff(bin_stat) / bin_stat[:-1]) * 100])
    keepgoing = pctdiff > 10  # keep going on 10% differences

    while any(keepgoing):
        # new bin edeges, splitting up any with > tol% difference
        new_bin_edges = []
        for i, v in enumerate(pctdiff):
            if v < pcttol:
                new_bin_edges.append(bin_edges[i])
            elif (
                len(binnumber[binnumber == i]) < minpnts
            ):  # too few points to split
                new_bin_edges.append(bin_edges[i])
            else:
                new_bin_edges.append((bin_edges[i] + bin_edges[i - 1]) / 2)
                new_bin_edges.append(bin_edges[i])
        new_bin_edges.append(bin_edges[-1])

        # recalculate bins
        _stats = binned_stats(x, y, statistic="median", bins=new_bin_edges)
        bin_stat, bin_edges, binnumber = _stats
        pctdiff = np.array(
            [
                0,
                *np.abs(np.diff(bin_stat) / bin_stat[:-1]) * 100,
            ]  # (0% diff with self)
        )

        # test if should keep re-binning
        keepgoing = []
        for i, v in enumerate(pctdiff):
            if v < 10:
                keepgoing.append(False)
            elif len(binnumber[binnumber == i]) < 25:
                keepgoing.append(False)
            else:
                keepgoing.append(True)
    # Bin Centers
    bin_cents = bin_edges[:-1] + np.diff(bin_edges) / 2

    # Adding back units
    if xunits is not None:
        bin_cents = bin_cents * xunits
        bin_edges = bin_edges * xunits
    if yunits is not None:
        bin_stat = bin_stat * yunits

    return bin_stat, bin_cents, bin_edges, binnumber, pctdiff


# /def


#############################################################################
# END
