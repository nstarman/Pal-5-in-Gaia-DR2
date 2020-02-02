#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : initializing selection functions
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""logging initialization."""

__author__ = "Nathaniel Starkman"

##############################################################################
# IMPORTS

# GENERAL
from scipy.interpolate import interp1d
from astropy import units as u

# CUSTOM
from astroPHD import LogFile

# PROJECT-SPECIFIC
from .select import inRange


##############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


##############################################################################
# CODE
##############################################################################


def select_shift_CMD(
    data,
    isochrone,
    lhs=0.0,
    rhs=0.08,
    low=0.035,
    up=0.0385,
    isorng=[18, 24],
    datarng=[18.1, 23.5],
    g_name="g dx",
    gmr_name="g-r dx",
    iso_g_name="g PS",
    iso_gmr_name="g-r PS",
    fill_value="extrapolate",
    logger=_LOGFILE,
    verbose=None,
):
    """Shifted CMD selection.

    Parameters
    ----------
    data : QTable
        the sky datatable
    isochrone : QTable
        isochrone datatable
    lhs : scalar
        shift the isochrone left
    rhs : scalar
        shift the isochrone right
        rhs > lhs
    low : scalar
        shift the isochrone down
    up : scalar
        shift the isochrone up
    isorng : (2,) list
        restrict the range applied on the isochrone
    datarng : (2,) list
        restrict the range applied on the data
    g_name : str, optional  (default 'g dx')
        name of the g column in *data*
    gmr_name : str, optional  (default 'g-r dx')
        name of the g-r column in *data*
    iso_g_name : str, optional  (default 'g PS')
        name of the g column in *isochrone*
    iso_gmr_name : str, optional  (default 'g-r PS')
        name of the g-r column in *isochrone*
    fill_value : str, optional  (default 'extrapolate')
        interp1d fill_value for lhs & rhs

    # Logging
    logfile : logger, optional
    verbose : int, optional
        the degree of verbosity
        None) (default): use instantiated value
        0) None; 1) status report, >=2) step-by-step

    """
    # -----------------------------------------------------
    # Checking

    assert lhs < rhs
    assert low < up

    # -----------------------------------------------------

    # restricting isochrone range
    _iso_ind = inRange(isochrone[iso_g_name], rng=isorng * u.mag)

    # isochrone lhs, and left shift
    isocr_spl_lhs = interp1d(
        isochrone[iso_g_name][_iso_ind] + lhs * u.mag,
        isochrone[iso_gmr_name][_iso_ind],
        fill_value=fill_value,
    )

    # isochrone rhs, and right shift
    isocr_spl_rhs = interp1d(
        isochrone[iso_g_name][_iso_ind] + rhs * u.mag,
        isochrone[iso_gmr_name][_iso_ind],
        fill_value=fill_value,
    )

    # -----------------------------------------------------

    # restricting data range
    _data_ind = inRange(data[g_name], rng=datarng * u.mag)

    # evaluating isochrone splines
    # shifting low and up
    evl_spl_MSp_low = isocr_spl_lhs(data[g_name][_data_ind]) + low
    evl_spl_MSp_up = isocr_spl_rhs(data[g_name][_data_ind]) + up
    rng = [evl_spl_MSp_low, evl_spl_MSp_up] * u.mag

    # -----------------------------------------------------

    # CMD
    shiftCMD = inRange(data[g_name], rng=datarng * u.mag)
    _cmd_ind = inRange(data[gmr_name][_data_ind], rng=rng)
    shiftCMD[shiftCMD] = _cmd_ind

    # -----------------------------------------------------
    # report

    logger.report(
        f"made select_shift_CMD",
        f"select_shift_CMD:\n\tlhs={lhs}, rhs={rhs}\n\tlow={low}, up={up}\
          \n\tisorng={isorng}, datarng={datarng}\
          \n\tg_name={g_name}, gmr_name={gmr_name}\
          \n\tiso_g_name={iso_g_name}, iso_gmr_name={iso_gmr_name}\
          \n\tfill_value={fill_value}",
        verbose=verbose,
    )

    # -----------------------------------------------------

    return shiftCMD


# /def

##############################################################################
# END
