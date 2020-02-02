#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : g selection
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""g selection."""

__author__ = "Nathaniel Starkman"


##############################################################################
# IMPORTS

# GENERAL
from astropy import units as u

# CUSTOM
from astroPHD import ObjDict, LogFile

# PROJECT-SPECIFIC
from .select import inRange

##############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


##############################################################################
# CODE
##############################################################################


def select_g_range(
    df, low=20, up=20.7, g_name="g dx", logger=_LOGFILE, verbose=None
):
    """select_g_range.

    Parameters
    ----------
    df : (Q)Table

    low : float
        units of magnitudes
    up : float
        units of magnitudes

    """

    res = inRange(df[g_name], rng=[low, up] * u.mag)

    # -----------------------------------------------------
    # report

    logger.report(
        "Made g range Selection",
        f"select_g_range:\n\tlhs={low}, up={up}\
          \n\tg_name={g_name}",
        verbose=verbose,
    )

    # -----------------------------------------------------

    return res


# /def

##############################################################################
# END
