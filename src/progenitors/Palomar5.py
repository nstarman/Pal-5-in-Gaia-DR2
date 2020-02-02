# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#
# TITLE   : Palomar5
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# -----------------------------------------------------------------------------

# Docstring and Metadata
"""Palomar 5 Information."""

__author__ = "Nathaniel Starkman"

##############################################################################
# IMPORTS

# Astropy
from astropy import units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord

# PROJECT-SPECIFIC
from astroPHD import ObjDict

##############################################################################
# PARAMETERS

_mas_yr = u.mas / u.yr


##############################################################################
# CODE
##############################################################################

progenitor = ObjDict(
    "Palomar 5",
    coords={
        "main": SkyCoord(
            ra=229.018 * u.deg,
            dec=-0.124 * u.deg,
            distance=23.2 * u.kpc,
            pm_ra_cosdec=-2.296 * _mas_yr,
            pm_dec=-2.257 * _mas_yr,
            radial_velocity=-58.7 * u.km / u.s,
        ),
        "vasiliev": SkyCoord(
            ra=229.022 * u.deg,
            dec=-0.223 * u.deg,
            distance=23.2 * u.kpc,
            pm_ra_cosdec=-2.736 * _mas_yr,
            pm_dec=-2.646 * _mas_yr,
            radial_velocity=-58.6 * u.km / u.s,
        ),
    },
    info=QTable(
        {
            "angular size ra": [13.18] * u.arcmin,
            "angular size dec": [11.48] * u.arcmin,
            "tidal radius": [0.145] * u.kpc,
            "tidal radius error": [0.009] * u.kpc,
            "tidal radius angular": [21.2] * u.arcmin,
            "tidal radius angular error": [None],
            "half-mass radius": [2.73] * u.arcmin,
            "pmra": [-2.296] * _mas_yr,
            "pmra_err": [0.3] * _mas_yr,
            "pmdec": [-2.257] * _mas_yr,
            "pmdec_err": [0.3] * _mas_yr,
        }
    ),
    references=[
        "main coord from Bovy prog Orbit",
        "angular size from SIMBAD",
        "tidal radius from Ibata 2017",
        "half-mass radius from http://physwww.mcmaster.ca/~harris/mwgc.dat",
        "proper motion from bovy orbit",
    ],
)

progenitor["t+"] = {"end": 0.1 * u.Gyr, "num": 10000}
progenitor["t-"] = {"end": -0.1 * u.Gyr, "num": 10000}
