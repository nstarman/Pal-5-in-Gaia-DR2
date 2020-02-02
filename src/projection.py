# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : projection
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Projection.

Routine Listings
----------------
get_data_along_orbit
select_stars_in_an_arm
digitize_star_along_stream

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "get_data_along_orbit",
    "select_stars_in_an_arm",
    "digitize_star_along_stream",
]


##############################################################################
# IMPORTS

# GENERAL
import numpy as np
from scipy.stats import binned_statistic

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable

# CUSTOM
from astroPHD import ObjDict, LogFile
from astroPHD.plot import starkplot as plt

# PROJECT-SPECIFIC
from .plot import plot_data_along_orbit, plot_data_along_orbit_in_window


##############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile

Table_like = (Table, QTable)


##############################################################################
# CODE
##############################################################################


def get_data_along_orbit(
    data,
    orbit,
    plot=True,
    frame="cust",
    lon_angle="phi1",
    lat_angle="phi2",
    **kw
):
    """get_data_along_orbit.

    Parameters
    ----------
    data : astropy SkyCoord
        the data
    orbit : astropy SkyCoord
    plot : bool

    Other Parameters
    ----------------
    kw
        passed to starkplot's plt.set
        see documentation for options.

    Returns
    -------
    idxres : ObjDict
        contains idxcatalog, idxo, arc, wsky

    """
    # match to catalog for distance reference
    idx2, d2d, _ = data.match_to_catalog_sky(orbit, nthneighbor=1)

    # data subset
    sep_ind = d2d <= 3 * u.deg  # points within width of orbit
    idxo = idx2[sep_ind]  # index into orbit
    wsky = d2d[sep_ind]  # sky separation
    idxcatalog = np.where(sep_ind == np.True_)[0]  # index into data

    # getting distance along the arc of the orbit
    orbphi1 = getattr(orbit, lon_angle)
    orbphi2 = getattr(orbit, lat_angle)
    arc = np.zeros_like(orbphi1)
    arc[1:] = np.cumsum(
        np.sqrt(
            np.diff(orbphi1) ** 2 * np.cos(orbphi2[:-1]) ** 2
            + np.diff(orbphi2) ** 2
        )
    )

    # storing
    idxres = ObjDict(
        "Data Along Orbit",
        idxcatalog=idxcatalog,
        idxo=idxo,
        # along & perpendicular to arc
        arc=arc,
        wsky=wsky,
    )

    # plotting
    if plot:
        # 1st Plot: Closeup
        plot_data_along_orbit(
            data,
            orbit,
            idxres,
            frame=frame,
            lon_angle=lon_angle,
            lat_angle=lat_angle,
            **kw
        )

        # 2nd Plot: SkyWindow
        plot_data_along_orbit_in_window(
            data,
            orbit,
            idxres,
            frame=frame,
            lon_angle=lon_angle,
            lat_angle=lat_angle,
        )
        plt.plot(orbit.phi1, orbit.phi2)
    # /if

    return idxres


# /def


# -----------------------------------------------------------------------------


def select_stars_in_an_arm(
    data,
    orbit,
    idxres,
    skycut=Ellipsis,
    bins=61,
    num_side=1,
    digitize=True,
    numstars=20,
):
    """plot_data_orbit_projection.

    Parameters
    ----------
    data : SkyCoord, (Q)Table
    orbit : SkyCoord, (Q)Table
    idxres :
        returned by get_data_along_orbit
    skycut : index array, optional  (default Ellipsis)
    bins : number of bins
        must be odd

    Returns
    -------
    res: ObjDict
        .statistic
        .bin_meds
        .binnumber
        .cent_inds
        .instream_catalog
        .instream

    """
    # --------------------
    if not bool(bins % 2):  # it's even
        raise ValueError("bins must be odd")

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

    x = idxres.arc[idxo][skycut].to_value("deg")

    y = idxres.wsky.copy()
    y[dataphi2 < orbitphi2] *= -1
    y = y[skycut].to_value("deg")

    # --------------------

    stats, bin_edges, binnum = binned_statistic(
        y, x, bins=bins, statistic="count"
    )

    # medians of bins
    bin_meds = bin_edges[:-1] + np.diff(bin_edges) / 2

    # indices of the stream
    cent_ind = int(np.around(np.max(binnum) / 2))
    cent_inds = np.arange(num_side * 2 + 1) - num_side + cent_ind - 1

    # --------------------
    # backing out the indices of stars in the stream

    # inside skycut, inside catalog
    streaminskycut = np.zeros_like(binnum, dtype=bool)
    for ind in cent_inds:
        streaminskycut |= binnum == ind

    # inside catalog
    streamincatalog = skycut.copy()
    streamincatalog[streamincatalog] = streaminskycut

    # in the whole dataset
    instream = idxcatalog[streamincatalog]

    # --------------------

    res = ObjDict(
        "instream",
        statistic=stats,
        bin_meds=bin_meds,
        binnumber=binnum,
        cent_inds=cent_inds,
        instream_catalog=streamincatalog,
        instream=instream,
    )

    return res


# /def


# -----------------------------------------------------------------------------


def digitize_star_along_stream(data, instream, bins=10):
    """digitize_star_along_stream.

    Parameters
    ----------
    data: (Q)Table
    instream: array_like
        index / bool array

    Returns
    -------
    df: (n, 3) ndarray
        columns ra, dec, dec_err

    Notes
    -----
    dec_err estimated as bin_width / sqrt(numpoints)
    works with output of select_stars_in_an_arm

    """
    if isinstance(data, SkyCoord):
        dataphi1 = data.phi1
        ra, dec = data.icrs.ra, data.icrs.dec
    elif isinstance(data, (Table, QTable)):
        dataphi1 = data["phi1"]
        ra, dec = data["ra"], data["dec"]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    x = dataphi1[instream]

    bin_edges = np.histogram_bin_edges(x, bins=bins)  # binning
    binnums = np.digitize(x, bin_edges[:-1])  # getting binnumber

    print(len(x), len(binnums))
    print(bins, len(np.unique(binnums)))

    avg_ra = np.full(bins, np.nan)
    avg_dec = np.full(bins, np.nan)
    avg_dec_err = np.full(bins, np.nan)

    for i, b in enumerate(np.unique(binnums)):
        ind = binnums == b
        avg_ra[i] = ra[instream][ind].mean().value
        avg_dec[i] = dec[instream][ind].mean().value
        avg_dec_err[i] = np.diff(bin_edges)[i] / np.sqrt(
            ind.sum()
        )  # width/sqrt(numpoints)

    return np.c_[avg_ra, avg_dec, avg_dec_err]


# /def


##############################################################################
# END
