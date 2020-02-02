# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : mwdust_util
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""utilities for mwdust.

TODO use _LOGFILE inside load_dust_color
TODO use _LOGFILE inside load_dust_gri

"""

__author__ = "Nathaniel Starkman"

#############################################################################
# IMPORTS

# GENERAL
import sys
import numpy as np
import mwdust

# Astropy
from astropy import units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord

# CUSTOM
from astroPHD import LogFile


#############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


#############################################################################
# CODE
#############################################################################


def load_dust_color(
    color,
    fname,
    df=None,
    distance=None,
    save=True,
    logger=_LOGFILE,
    verbose=None,
):
    """PanSTARRS1 dust color.

    Parameters
    ----------
    color
    fname
    df
    distance
    save
    logger
    verbose

    Returns
    -------
    dust_c

    """
    try:
        dust_c = np.load(fname.format(color)) * u.mag

    except FileNotFoundError:
        if (df is None) or (distance is None):
            raise ValueError(
                f"load_dust_color({color}) needs *df* & *distance*"
            )

        sys.stdout.write(f"could not load {color}")
        sys.stdout.flush()

        SFDc = mwdust.SFD(filter=f"PS1 {color}")
        dust_c = SFDc(df["l"], df["b"], distance.to_value(u.kpc)) * u.mag

        if save:
            dust_c.value.dump(fname.format(color))

        sys.stdout.write(f"\033 -> finised {color}, ")
        sys.stdout.flush()

    else:
        sys.stdout.write(f"\033 {color}, ")
        sys.stdout.flush()

    return dust_c


# /def


# ----------------------------------------------------------------------------


def load_dust_gri(
    fname,
    df=None,
    distance=None,
    recalculate=False,
    save=True,
    save_intermediate=False,
    logger=_LOGFILE,
    verbose=None,
):
    """Load_dust_gri.

    panstarrs1 dust in g,r,i
    tries to load from fname.format('gri'), else tries to make table

    Arguments
    ---------
    fname:
        must support color substitution
        ex: f'../data/nbody3/skydata/ps1dust_040112_{}.dat'
    df: QTable
        table of coordinates.
        only used if load fails
    distance: quantity
        the distance at which to calculate the dust extinction
        only used if load fails
    save: bool
        whether to save the table
        only used if load fails
    save_intermediate:
        whether to save the intermediate tables
        only used if load fails

    Returns
    -------
    ps1dust_gri: QTable
        table of dust extinction
        [l, b], g, r, i

    """
    df = QTable(df)  # just making sure

    try:
        if recalculate:
            raise FileNotFoundError

        ps1dust_gri = QTable.read(fname.format("gri"), format="ascii.ecsv")

    except FileNotFoundError:
        sys.stdout.write("\ncould not load gri table")
        sys.stdout.flush()

        # loading g
        ps1dust_g = load_dust_color(
            "g",
            fname,
            df=df,
            distance=distance,
            save=save_intermediate,
            logger=logger,
            verbose=verbose,
        )

        ps1dust_r = load_dust_color(
            "r",
            fname,
            df=df,
            distance=distance,
            save=save_intermediate,
            logger=logger,
            verbose=verbose,
        )

        ps1dust_i = load_dust_color(
            "i",
            fname,
            df=df,
            distance=distance,
            save=save_intermediate,
            logger=logger,
            verbose=verbose,
        )

        coord = SkyCoord(
            l=df["l"], b=df["b"], distance=distance, frame="galactic"
        )

        ps1dust_gri = QTable(
            [coord, ps1dust_g, ps1dust_r, ps1dust_i],
            names=("coord", "g", "r", "i"),
        )

        if save:
            ps1dust_gri.write(fname.format("gri"), format="ascii.ecsv")

        sys.stdout.write("\033" + " -> assembled gri table, ")
        sys.stdout.flush()

    else:
        # check if need to recalculate the dust because doesn't match df
        if len(df["g"]) != len(ps1dust_gri["g"]):
            logger.write("need to recalculate dust")
            ps1dust_gri = load_dust_gri(
                fname,
                df=df,
                distance=distance,
                recalculate=True,
                save=save,
                save_intermediate=save_intermediate,
                logger=logger,
                verbose=verbose,
            )

        sys.stdout.write("loaded gri table")
        sys.stdout.flush()

    return ps1dust_gri


# /def


##############################################################################
# END
