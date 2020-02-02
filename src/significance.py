# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : compute significance
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Compute Stream Significance.

Routine Listings
----------------
calculate_significance


"""

__author__ = "Nathaniel Starkman"

__all__ = ["calculate_significance"]


##############################################################################
# IMPORTS

# GENERAL
import numpy as np
from scipy.stats import poisson

from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable

import matplotlib.pyplot as plt

from typing import Optional, Union

# CUSTOM
from astroPHD import ObjDict


##############################################################################
# PARAMETERS

data_typing = Union[SkyCoord, Table]
array_like = Union[list, tuple, np.ndarray]


##############################################################################
# CODE
##############################################################################


def calculate_significance(
    data: data_typing,
    orbit: data_typing,
    idxres: Optional[ObjDict],
    skycut: array_like,
    nsample: int = 27,
    cdf: array_like = [],
    plot=True,
):
    """Calculate significance.

    Parameters
    ----------
    data: SkyCoord or (Q)Table
    orbit: SkyCoord or (Q)Table
    idxres: ObjDict, optional
    skycut: array_like
    nsample: int, optional
    cdf: array_lik, optional

    Returns
    -------
    significance: float
        the significance of the stream detection

    """
    # --------------------
    if idxres is None:  # assumed already used idxcatalog & idxo
        idxcatalog = np.ones(len(data), dtype=bool)
        idxo = np.ones(len(data), dtype=bool)
    else:
        idxcatalog = idxres.idxcatalog
        idxo = idxres.idxo

    # --------------------

    if isinstance(data, SkyCoord):
        # dataphi1 = data.phi1[idxcatalog]
        dataphi2 = data.phi2[idxcatalog]
    elif isinstance(data, (Table, QTable)):
        # dataphi1 = data['phi1'][idxcatalog]
        dataphi2 = data["phi2"][idxcatalog]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    if isinstance(orbit, SkyCoord):
        orbitphi2 = orbit.phi2[idxo]
    elif isinstance(orbit, (Table, QTable)):
        orbitphi2 = orbit["phi2"][idxo]
    else:
        raise TypeError("orbit is not a SkyCoord or (Q)Table")

    # --------------------

    x = idxres.arc[idxres.idxo][skycut].to_value("deg")

    y = idxres.wsky.copy()
    y[dataphi2 < orbitphi2] *= -1
    y = y[skycut].to_value("deg")

    # --------------------

    if plot:

        # plot 1
        fig, ax = plt.subplots(1, 2, figsize=(6, 3))
        ax[0].scatter(x, y)
        freq, ybin, _ = ax[1].hist(y, bins=nsample)
        ax[1].axhline(np.mean(freq), c="k", ls=":")

        # plot 2
        mu = np.median(freq)
        sample = poisson.rvs(mu, size=nsample)

        plt.figure()
        plt.hist(y, bins=nsample)
        plt.axhline(np.mean(freq), c="k", ls=":")
        plt.scatter(ybin[:-1], sample, c="k", zorder=2)

    # /if

    if not cdf:
        raise Exception
    else:
        significance = (nsample - 2) * np.prod(
            [(1 - poisson.cdf(x, mu)) for x in cdf]
        )

    print(f"p value = {significance * 100:.10f} %")

    return significance


# /def

##############################################################################
# END
