#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   :
# AUTHOR  : Nathaniel Starkman
# PROJECT :
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""combine selections."""

__author__ = "Nathaniel Starkman"


###############################################################################
# IMPORTS

# GENERAL
from astropy import units as u

# CUSTOM
from astroPHD import LogFile

# PROJECT-SPECIFIC
from .cmd import select_shift_CMD
from .pm import select_pm_circle
from .gsel import select_g_range


###############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False, verbose=0)


###############################################################################
# CODE
##############################################################################


def select_pm_cmd_g_cut(
    data,
    isochrone,
    # PM
    pm_x0_lon=-2.5,
    pm_x0_lat=-2.5,
    pm_dx_lon=3,
    pm_dx_lat=3,
    # CMD
    cmd_lhs=0.0,
    cmd_rhs=0.08,
    cmd_low=0.035,
    cmd_up=0.0385,
    # G
    g_low=20,
    g_up=20.7,
    # kwargs for PM, CMD, G
    lon_name="pmra",
    lat_name="pmdec",  # pm
    lon_units=u.mas / u.yr,
    lat_units=u.mas / u.yr,  # pm
    g_name="g dx",
    gmr_name="g-r dx",  # cmd, g
    iso_g_name="g PS",
    iso_gmr_name="g-r PS",  # cmd
    isorng=[18, 24],
    datarng=[18.1, 23.5],  # cmd
    fill_value="extrapolate",  # cmd
    # logging
    logger=_LOGFILE,
    verbose=None,
    **kw
):
    """select_pm_cmd_g_cut.

    Parameters
    ----------
    df : QTable
        the sky datatable
        must contain columns *pm_lon*, *pm_lat* that are in units
        *lon_units*, *lat_units*, respectively.
    isochrone : QTable
        isochrone datatable

    # PM
    pm_x0_lon : float
        proper motion xo in the longitude direction
        units are in *lon_units*
    pm_x0_lat : float
        proper motion xo in the latitude direction
        units are in *lat_units*
    pm_dx_lon : float
        proper motion dx in the longitude direction
        units are in *lon_units*
    pm_dx_lat : float
        proper motion dx in the latitude direction
        units are in *lat_units*

    # CMD
    cmd_lhs : scalar
        shift the isochrone left
    cmd_rhs : scalar
        shift the isochrone right
        rhs > lhs
    cmd_low : scalar
        shift the isochrone down
    cmd_up : scalar
        shift the isochrone up

    # kwargs for PM, CMD, G
    lon_name : str, optional  (default 'pmra')
        name of the longitudinal direction
    lat_name : str, optional  (default 'pmdec')
        name of the latitudinal direction
    lon_units : astropy units, optional  (default u.mas / uyr)
        units of the longitudinal direction
    lat_units : astropy units, optional  (default u.mas / uyr)
        units of the latitudinal direction
    g_name : str, optional  (default 'g dx')
        name of the g column in *data*
    gmr_name : str, optional  (default 'g-r dx')
        name of the g-r column in *data*
    iso_g_name : str, optional  (default 'g PS')
        name of the g column in *isochrone*
    iso_gmr_name : str, optional  (default 'g-r PS')
        name of the g-r column in *isochrone*
    isorng : (2,) list
        restrict the range applied on the isochrone
    datarng : (2,) list
        restrict the range applied on the data
    fill_value : str, optional  (default 'extrapolate')
        interp1d fill_value for lhs & rhs

    # Logging
    logfile : logger, optional
    verbose : int, optional
        the degree of verbosity
        None) (default): use instantiated value
        0) None; 1) status report, >=2) step-by-step

    """
    # showing unused kwargs for safety (in case write in wrong kw)
    # have this to absorb lmfit extra parameters
    logger.report("unused kwargs: {}".format(kw), verbose=verbose)

    # proper motion
    selpm = select_pm_circle(
        data,
        x0_lon=pm_x0_lon,
        x0_lat=pm_x0_lat,
        dx_lon=pm_dx_lon,
        dx_lat=pm_dx_lat,
        lon_name=lon_name,
        lat_name=lat_name,
        lon_units=lon_units,
        lat_units=lat_units,
        logger=_LOGFILE,
        verbose=None,
    )

    # cmd
    selcmd = select_shift_CMD(
        data,
        isochrone,
        lhs=cmd_lhs,
        rhs=cmd_rhs,
        low=cmd_low,
        up=cmd_up,
        isorng=isorng,
        datarng=datarng,
        g_name=g_name,
        gmr_name=gmr_name,
        iso_g_name=iso_g_name,
        iso_gmr_name=iso_gmr_name,
        fill_value=fill_value,
        logger=_LOGFILE,
        verbose=None,
    )

    # g
    selg = select_g_range(
        data, low=g_low, up=g_up, g_name=g_name, logger=_LOGFILE, verbose=None
    )

    # combining
    comb_sel = selpm & selcmd & selg

    return comb_sel


# /def

###############################################################################
# END
