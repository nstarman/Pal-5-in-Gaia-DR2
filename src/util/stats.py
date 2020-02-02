# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#
# TITLE   : coordinates initialization file
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# -----------------------------------------------------------------------------

# Docstring and Metadata
"""Statistics on SkyCoords."""

__author__ = "Nathaniel Starkman"

__all__ = [
    "mean"
]


#############################################################################
# IMPORTS

# GENERAL
from astropy.coordinates import SkyCoord


#############################################################################
# CODE
#############################################################################


def mean(sc):
    """Average of skycoords.

    Parameters
    ----------
    sc : SkyCoord
        non-scalar SkyCoord

    Returns
    -------
    avg_sc : SkyCoord
        average of `sc`

    """
    frame = sc.frame
    representation_type = sc.representation_type

    x, y, z = sc.icrs.cartesian.xyz.mean(axis=1)
    v_x, v_y, v_z = sc.icrs.velocity.d_xyz.mean(axis=1)

    avg_sc = SkyCoord(
        x=x,
        y=y,
        z=z,
        v_x=v_x,
        v_y=v_y,
        v_z=v_z,
        frame="icrs",
        representation_type="cartesian",
    )
    avg_sc.representation_type = representation_type

    return avg_sc.transform_to(frame)


# /def

#############################################################################
# END
