# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : isochrone
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""isochrone."""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# GENERAL
import numpy as np
from astropy import units as u
from astropy.table import Table


# PROJECT-SPECIFIC
from .units_decorators import quantity_io
from . import table_utils as mqt


#############################################################################
# CODE
#############################################################################

#############################################################################
# Distance Modulus


@quantity_io(d="length", A=u.mag, annot2dfu=True, default_units={"A": u.mag})
def distanceModulus_magnitude(d: u.pc, A=0 * u.mag, obs=True, **kw) -> u.mag:
    """Distance Modulus.

    equation:  DM = 5 log10(d / 10) + A

    mtrue - M = 5 log10(d / 10)
    if there is line-of-sight extinction
    mobs = mtrue + A
    | mobs - M = 5 log10(d / 10) + A
    | true - M = 5 log10(d / 10) - A

    Arguments
    ---------
    d: scalar, array, Quantity
        distance
        no units => parsecs
    A: scalar, array, Quantity (in mag)
        extinction in magnitudes
    obs: bool (default True)
        whether to return (mobs - M) or (mtrue - M)
        defualt: (mobs - M)
        **don't change unless specifically required

    Returns
    -------
    DM: scalar, array
        default units: u.mag

    """
    if not obs:
        A *= -1

    return (5 * u.mag) * (np.log10(d.to_value(u.pc)) - 1) + A


# /def


#############################################################################
# Table Readers


def CMD3p1TableRead(
    fname, detector="ps1", fmat="ascii", header_start=7, distance=None
):
    """CMD3p1TableRead.

    Read a CMD table from  CMD 3.1 @ http://stev.oapd.inaf.it/cgi-bin/cmd
    CMD starts in absolute magnitudes

    """
    isocr = Table.read(fname, format=fmat, header_start=header_start)
    isocr.meta["comments"] = isocr.meta["comments"][: header_start - 1]

    if detector in ("cfht",):
        absmagcols = ("Umag", "Gmag", "Rmag", "Imag", "Zmag")
        appmagcols = ("u", "g", "r", "i", "z")
        mqt.rename_columns(
            isocr,
            ("u*mag", "Umag"),
            ("g'mag", "Gmag"),
            ("r'mag", "Rmag"),
            ("i'mag", "Imag"),
            ("z'mag", "Zmag"),
        )

    elif detector in ("ps1", "panstarrs", "panstarrs1"):
        absmagcols = ("Gmag", "Rmag", "Imag", "Zmag", "Ymag", "Wmag")
        appmagcols = ("g", "r", "i", "z", "y", "w")
        mqt.rename_columns(
            isocr,
            ("gP1mag", "Gmag"),
            ("rP1mag", "Rmag"),
            ("iP1mag", "Imag"),
            ("zP1mag", "Zmag"),
            ("yP1mag", "Ymag"),
            ("wP1mag", "Wmag"),
        )

    udict = {
        "Age": u.yr,
        "Mass": u.Msun,
        "mbolmag": u.mag,
        "Umag": u.mag,
        "Gmag": u.mag,
        "Rmag": u.mag,
        "Imag": u.mag,
        "Zmag": u.mag,
        "Ymag": u.mag,
        "Wmag": u.mag,
    }

    isocr = mqt.add_units_to_Table(isocr, udict=udict)  # making QTable

    if set(["Gmag", "Rmag"]).issubset(isocr.colnames):
        mqt.add_color_col(isocr, "Gmag", "Rmag", color="g-r")

    # Adjusting to apparent magnitude
    if distance is not None:
        DM = distanceModulus_magnitude(distance)
        appmags = [isocr[n] + DM for n in absmagcols]

        isocr.add_columns(appmags, names=appmagcols)
    else:
        print("Not adding apparent mags (cols ugriz)")

    return isocr

# /def


# -------------------------------------------------------------------------


def readCMDTablestevoapd(
    fname, detector="ps1", fmat="ascii", header_start=7, distance=None
):
    """Calls CMD3p1TableRead."""
    return CMD3p1TableRead(
        fname,
        detector=detector,
        fmat=fmat,
        header_start=header_start,
        distance=distance,
    )

# /def


##############################################################################
# END
