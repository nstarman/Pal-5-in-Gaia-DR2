# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : util initialization
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""util.

Routine Listings
----------------
astrarray

"""

__author__ = "Nathaniel Starkman"

#############################################################################
# IMPORTS

# GENERAL
import numpy as np
from collections import OrderedDict

# CUSTOM
from astroPHD import LogFile, ObjDict

# PROJECT-SPECIFIC
from .pickle import dump as _dump, load as _load


#############################################################################
# CODE
#############################################################################


def astrarray(arr):
    """Quantity array.

    converts a list of quantities to a Quantity list

    Parameters
    ----------
    arr: array_like

    Examples
    --------
    >>> astrarray([0*u.deg, 1*u.deg])
    [0, 1] * u.deg

    """
    astrclass = arr[0].__class__
    unit = arr[0].unit

    return astrclass(np.asarray([a.to_value(unit) for a in arr]) * unit)


# /def


#############################################################################
# END
