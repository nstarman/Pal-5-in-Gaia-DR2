# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : initializing selection functions
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Proper motion selections."""

__author__ = "Nathaniel Starkman"

##############################################################################
# IMPORTS

# GENERAL
from astropy import units as u

# CUSTOM
from astroPHD import LogFile

# Project-SPECIFIC
from .select import ellipse


##############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


##############################################################################
# CODE
##############################################################################


##############################################################################
# Select Proper Motion


def select_pm_circle(
    df,
    x0_lon=-2.5,
    x0_lat=-2.5,
    dx_lon=3,
    dx_lat=3,
    lon_name="pmra",
    lat_name="pmdec",
    lon_units=u.mas / u.yr,
    lat_units=u.mas / u.yr,
    logger=_LOGFILE,
    verbose=None,
):
    """Select proper motion circle.

    Parameters
    ----------
    df : QTable
        the sky datatable
        must contain columns *pm_lon*, *pm_lat* that are in units
        *lon_units*, *lat_units*, respectively.

    x0_lon : float
        xo in the longitude direction
        units are in *lon_units*
    x0_lat : float
        xo in the latitude direction
        units are in *lat_units*
    dx_lon : float
        dx in the longitude direction
        units are in *lon_units*
    dx_lat : float
        dx in the latitude direction
        units are in *lat_units*
    lon_name : str, optional  (default 'pmra')
        name of the longitudinal direction
    lat_name : str, optional  (default 'pmdec')
        name of the latitudinal direction
    lon_units : astropy units, optional  (default u.mas / uyr)
        units of the longitudinal direction
    lat_units : astropy units, optional  (default u.mas / uyr)
        units of the latitudinal direction

    # Logging
    logfile : logger, optional
    verbose : int, optional
        the degree of verbosity
        None) (default): use instantiated value
        0) None; 1) status report, >=2) step-by-step

    """
    pmcirc = ellipse(
        df[lon_name].to_value(lon_units),
        df[lat_name].to_value(lat_units),
        x0=(x0_lon, x0_lat),
        dx=(dx_lon, dx_lat),
    )

    # -----------------------------------------------------
    # report

    logger.report(
        "Made PM Selection",
        f"select_pm_circle:\n\tx0_lon={x0_lon}, x0_lat={x0_lat}\
          \n\tdx_lon={dx_lon}, dx_lat={dx_lat}\
          \n\tlon_name={lon_name}, lat_name={lat_name}\
          \n\tlon_units={lon_units}, lat_units={lat_units}",
        verbose=verbose,
    )

    # -----------------------------------------------------

    return pmcirc


# /def

##############################################################################
# END
