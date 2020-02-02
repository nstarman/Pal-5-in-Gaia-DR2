# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : plotfuncs
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""plot functions.

Routine Listings
----------------
plot_sky_window
plot_proper_motions
plot_data_along_orbit
plot_data_along_orbit_in_window
plot_data_orbit_projection
plot_radec_track_residual

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "plot_sky_window",
    "plot_proper_motions",
    "plot_data_along_orbit",
    "plot_data_along_orbit_in_window",
    "plot_data_orbit_projection",
    "plot_radec_track_residual",
]

##############################################################################
# IMPORTS

# GENERAL
import numpy as np

from scipy.ndimage.filters import gaussian_filter

from matplotlib import colors
from matplotlib.patches import Ellipse

from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable

from typing import Optional

# CUSTOM
from astroPHD.plot import starkplot as plt, mpl_decorator

# PROJECT-SPECIFIC
from .select import inRange


##############################################################################
# CODE
##############################################################################


@mpl_decorator(
    fig="new",
    figsize=(10, 5),
    xlabel=r"$\phi_1$ (deg)",
    ylabel=r"$\phi_2$ (deg)",
    aspect="equal",
    invert_axis="xyz",
    xkw={"fontsize": 15, "title": {"fontsize": 20}},
)
def plot_sky_window(
    data,
    ind=...,
    step: int = 0.05,
    sigma: tuple = (2, 2),
    nbody: Optional[dict] = None,
    indwinstr=None,
    orb: Optional[SkyCoord] = None,
    orbrng=None,
    frame="cust",
    lon_angle="phi1",
    lat_angle="phi2",
    cmap="Spectral",
    **kw,
):
    """Plot the sky data window.

    Parameters
    ----------
    data: astropy QTable
        the dataframe containing the sky window

    ind: Ellipse or array_like, optional
    step: float or (2,) list of floats, optional
        (default = 0.05)
        the number of degrees by which to bin the sky data
        if (2,) list of floats then taken as (xstep, ystep)
    sigma: float or (2,) list of floats, optional
        (default = 2)
        the sigma for the Gaussian filter
        if (2,) list of floats then taken as (xstep, ystep)

    nbody: dictionary, optional
        dictionary must have 'cust', which is a SkyCoord
    orb: SkyCoord, optional
    orbrng: array_like, optional

    frame: str, optional
    lon_angle: str, optional
    lat_angle: str, optional

    cmap:  str or matplotlib colormap
        the colormap to use for the Gaussia-filtered image
        TODO some reason need to specify, can't be in mpl_decorator?

    nbody_cmap:  (default = 'nipy_spectral_r')

    TODO
    ----
    explain normed=True, interpolation='nearest'
    support not only phi1, phi2

    Returns
    -------
    im: pyplot imshow
        Gaussian-filtered image

    """
    plt.grid(False)

    # --------------------
    # step
    if np.isscalar(step):
        xstep = ystep = step
    else:
        xstep, ystep = step

    # sigma
    if np.isscalar(sigma):
        xsigma = ysigma = sigma
    else:
        xsigma, ysigma = sigma

    # --------------------

    # getting x, y data
    if isinstance(data, SkyCoord):
        x = getattr(data, lon_angle)[ind].to_value("deg")
        y = getattr(data, lat_angle)[ind].to_value("deg")
    elif isinstance(data, (Table, QTable)):
        x = data[lon_angle][ind].to_value("deg")
        y = data[lat_angle][ind].to_value("deg")
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    # bins
    xbins = np.arange(x.min(), x.max() + xstep, xstep)
    ybins = np.arange(y.min(), y.max() + ystep, ystep)

    # making histogram
    # normalizing the counts
    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins], normed=False)
    H = np.rot90(H)  # rotating to correct orientation
    # getting the figure edges
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    # gaussian filter
    fH = gaussian_filter(H, sigma=sigma, order=0)

    # plotting
    im = plt.imshow(
        fH,
        extent=extent,
        interpolation="nearest",
        cmap=cmap,
        norm=colors.LogNorm(),
    )

    # -----------------------------------------------------------
    # nbody
    if nbody is not None:
        plt.scatter(
            getattr(nbody[frame], lon_angle)[indwinstr][::100],
            getattr(nbody[frame], lat_angle)[indwinstr][::100],
            c=nbody[frame].distance[indwinstr][::100].distmod,
            s=1,
            colorbar=True,
            cmap=kw.get("nbody_cmap", "nipy_spectral_r"),
        )

    # -----------------------------------------------------------
    # progenitor orbit
    if orb is not None:
        if orbrng is None:
            raise ValueError("orbrng cannot be <None>")

        orbwin = inRange(
            getattr(orb[frame], lon_angle),
            getattr(orb[frame], lat_angle),
            rng=orbrng,
        )

        plt.plot(
            getattr(orb[frame], lon_angle)[orbwin],
            getattr(orb[frame], lat_angle)[orbwin],
        )

    return im


# /def


# -----------------------------------------------------------------------------


@mpl_decorator(
    fig="new",
    figsize=(10, 10),
    title=(r"Pal5 PM", dict(fontsize=24)),
    xlabel=r"$\mu_{\alpha^*} \ (\rm{mas} \, \rm{yr}^{-1})$",
    ylabel=r"$\mu_\delta \ (\rm{mas} \, \rm{yr}^{-1})$",
    unit_labels=False,
    xlim=(-15, 10),
    ylim=(-15, 10),
    tight_layout=True,
    legend=dict(fontsize=15),
    xkw=dict(fontsize=20),
)
def plot_proper_motions(
    data,
    idx,
    bakidx: slice = slice(None, None, 100),
    orbit: Optional[SkyCoord] = None,
    orbitinwindow: Optional[list] = None,
    cutspm: Optional[list] = None,
    **kw,
):
    """Plot proper motions.

    Parameters
    ----------
    data: SkyCoord, (Q)Table
        must have pmra & pmdec columns
        if SkyCoord: .pm_ra_cosdec (or .pm_ra) & .pm_dec
        if (Q)Table: ['pmra'], ['pmdec']

    idx : slicer
        index or bool array to select the rows from `data`
    bakidx : slicer, optional
        background subsampling
    orbit : optional
        (default None)
    orbitinwindow : optional
        (default None)
        only used if orbit is not None
    cutspm : optional
        (default None)

    Returns
    -------
    line : scatterplot result
        plt.scatter(pmra, pmdec, **kw)

    """
    # --------------------

    if isinstance(data, SkyCoord):
        try:
            datapmra = data.pm_ra_cosdec
        except Exception:
            datapmra = data.pm_ra
        datapmdec = data.pm_dec
    elif isinstance(data, (Table, QTable)):
        datapmra = data["pmra"]
        datapmdec = data["pmdec"]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    notnan = ~(np.isnan(datapmra) | np.isnan(datapmdec))

    # Background
    _bakidx = np.zeros(len(datapmra), dtype=bool)
    _bakidx[bakidx] = True
    plt.scatter(datapmra[notnan & _bakidx], datapmdec[notnan & _bakidx], s=1)

    # Plotting
    _idx = np.zeros(len(datapmra), dtype=bool)
    _idx[idx] = True
    line = plt.scatter(datapmra[notnan & _idx], datapmdec[notnan & _idx], **kw)

    # Start value circle
    if cutspm is not None:
        ellipse = Ellipse(
            xy=(cutspm.x0_lon, cutspm.x0_lat),
            width=cutspm.dx_lon,
            height=cutspm.dx_lat,
            fill=False,
            edgecolor="k",
        )
        plt.gca().add_patch(ellipse)

    if orbit is not None:
        plt.plot(
            orbit["eq"].pm_ra_cosdec,
            orbit["eq"].pm_dec,
            label="all",
            lw=7,
            c="k",
        )

        if orbitinwindow is not None:
            plt.plot(
                orbit["eq"].pm_ra_cosdec[orbitinwindow],
                orbit["eq"].pm_dec[orbitinwindow],
                label="all",
                lw=7,
                c="tab:purple",
            )

        plt.scatter(
            orbit["orb"].pmra(t=None),
            orbit["orb"].pmdec(),
            label="prog",
            s=30,
            zorder=2,
        )

    return line


# /def


# -----------------------------------------------------------------------------


def plot_data_along_orbit(
    data, orbit, idxres, frame="cust", lon_angle="phi1", lat_angle="phi2", **kw
):
    """plot_data_along_orbit."""
    # --------------------

    if isinstance(data, SkyCoord):
        dataphi1 = getattr(data, lon_angle)[idxres.idxcatalog]
        dataphi2 = getattr(data, lat_angle)[idxres.idxcatalog]
    elif isinstance(data, (Table, QTable)):
        dataphi1 = data[lon_angle][idxres.idxcatalog]
        dataphi2 = data[lat_angle][idxres.idxcatalog]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    if isinstance(orbit, SkyCoord):
        orbitphi1 = getattr(orbit, lon_angle)[idxres.idxo]
        orbitphi2 = getattr(orbit, lat_angle)[idxres.idxo]
    elif isinstance(orbit, (Table, QTable)):
        orbitphi1 = orbit[lon_angle][idxres.idxo]
        orbitphi2 = orbit[lat_angle][idxres.idxo]
    else:
        raise TypeError("orbit is not a SkyCoord or (Q)Table")

    # --------------------

    # orbit
    plt.scatter(orbitphi1, orbitphi2, c="k", s=0.1)
    # data
    plt.scatter(dataphi1, dataphi2, s=0.1)
    # setting plot properties
    plt.set(
        aspect="equal",
        title=kw.pop("title", "Data Along Orbit"),
        xlabel=kw.pop("xlabel", rf"{lon_angle} [deg]"),
        ylabel=kw.pop("ylabel", rf"{lat_angle} [deg]"),
    )
    plt.set(**kw)  # any extra user keywords

    return


# /def


# -----------------------------------------------------------------------------


@mpl_decorator(fig="new")
def plot_data_along_orbit_in_window(
    data,
    orbit,
    idxres,
    windowcut=...,
    frame="cust",
    lon_angle="phi1",
    lat_angle="phi2",
    **kw,
):
    """plot_data_along_orbit_in_window.

    Parameters
    ----------
    data : SkyCoord
    orbit : SkyCoord

    """
    # --------------------

    if isinstance(data, SkyCoord):
        phi1 = getattr(data, lon_angle)[idxres.idxcatalog]
        phi2 = getattr(data, lat_angle)[idxres.idxcatalog]
    elif isinstance(data, (Table, QTable)):
        phi1 = data[lon_angle][idxres.idxcatalog]
        phi2 = data[lat_angle][idxres.idxcatalog]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    plot_sky_window(
        data,
        ind=...,
        orbrng=windowcut,
        fig=None,
        frame=frame,
        lon_angle=lon_angle,
        lat_angle=lat_angle,
    )
    plt.scatter(phi1, phi2, s=40, c="k", alpha=0.05)
    plt.set(**kw)

    return


# /def


# -----------------------------------------------------------------------------


def plot_data_orbit_projection(
    data,
    orbit,
    idxres,
    skycut=Ellipsis,
    filtered=True,
    step=0.07,
    sigma=(2, 2),
    add_sidehists=True,
    frame="cust",
    lon_angle="phi1",
    lat_angle="phi2",
    **kw,
):
    """plot_data_orbit_projection.

    Parameters
    ----------
    data : SkyCoord, (Q)Table
    orbit : SkyCoord, (Q)Table
    idxres :
        returned by get_data_along_orbit
    skycut : index array, optional  (default Ellipsis)
    filtered : bool, optional  (default True)
        whether to apply a gaussian filter
    step : float, optional  (default .07)
        the step size in the gaussian filter
    sigma : list, optional  (default (2, 2))
        the sigma in the guassian filter

    Returns
    -------
    im : matplotlib plot
        if not filtered:
            scatter(*, s=2, sidehists=True)
        if filtered:
            imshow(*, cmap='Spectral', norm=colors.LogNorm(),
                   interpolation='nearest')

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
        # dataphi1 = getattr(data, lon_angle)[idxcatalog]
        dataphi2 = getattr(data, lat_angle)[idxcatalog]
    elif isinstance(data, (Table, QTable)):
        # dataphi1 = data[lon_angle][idxcatalog]
        dataphi2 = data[lat_angle][idxcatalog]
    else:
        raise TypeError("data is not a SkyCoord or (Q)Table")

    # --------------------

    if isinstance(orbit, SkyCoord):
        # orbphi1 = getattr(orbit, lon_angle)
        orbphi2 = getattr(orbit, lat_angle)
    elif isinstance(orbit, (Table, QTable)):
        # orbphi1 = orbit[lon_angle]
        orbphi2 = orbit[lat_angle]
    else:
        raise TypeError("orbit is not a SkyCoord or (Q)Table")

    # --------------------

    # distance along the arc
    x = idxres.arc[idxo][skycut].to_value("deg")

    # perpendicular to arc
    y = idxres.wsky.copy()
    y[dataphi2 < orbphi2[idxo]] *= -1
    y = y[skycut].to_value("deg")

    # --------------------

    plt.grid(False)
    ax = plt.gca()

    if not filtered:
        im = plt.scatter(x, y, s=2, sidehists=add_sidehists)
    else:
        xbins = np.arange(x.min(), x.max() + step, step)
        ybins = np.arange(y.min(), y.max() + step, step)

        H, xedges, yedges = np.histogram2d(
            x, y, bins=[xbins, ybins], normed=True
        )
        H = np.rot90(H)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        fH = gaussian_filter(H, sigma=sigma, order=0)

        im = plt.imshow(
            fH,
            extent=extent,
            cmap="Spectral",
            norm=colors.LogNorm(),
            interpolation="nearest",
            xlabel=r"$s_1$",
            ylabel=r"$s_2$",
        )

        shargs = {k: v for k, v in kw.items() if k.startswith("sh")}
        [kw.pop(k) for k in shargs.keys()]
        plt.scatter(x, y, alpha=0, sidehists=add_sidehists, **shargs)

    axs = plt.gcf().axes
    axs[0].set_aspect("equal")
    axs[0].tick_params(
        axis="both", labelsize=13, color="black", colors="black"
    )
    axs[1].tick_params(axis="both", labelsize=12)
    axs[2].tick_params(axis="both", labelsize=12)
    axs[0].set_xlabel(r"$s_1 \ (\rm{deg})$", fontsize=20)
    axs[0].set_ylabel(r"$s_2 \ (\rm{deg})$", fontsize=20)
    if not isinstance(add_sidehists, bool):
        if "left" in add_sidehists:
            axs[0].yaxis.set_label_position("right")

    plt.set(figsize=(8, 3.5))
    plt.tight_layout()

    plt.set(**kw)
    plt.sca(ax)

    return im


# /def


# -----------------------------------------------------------------------------


def plot_radec_track_residual(data_radec, track_radec, track_interp, **kw):
    """plot_radec_track_residual.

    calculates residual as (dec - track_interp.dec) / dec_err

    Parameters
    ----------
    data_radec : (n, 3) array
        measured data.
        columns are [ra, dec, dec_err]
    track_radec :
    track_interp : function

    Returns
    -------
    fig : matplotlib figure
    (frame1, frame2) : matplotlib axes

    """
    # Preplotting
    ra, dec, dec_err = data_radec.T

    residual = dec - track_interp(ra)

    xmin, xmax = min(ra), max(ra)
    xlim = [xmin - 2.5, xmax + 2.5]

    ymin, ymax = min(dec), max(dec)
    ylim = [ymin - 2.5, ymax + 2.5]

    # Plotting
    fig = plt.figure()
    # main plot
    frame1 = fig.add_axes((0.1, 0.3, 0.8, 0.6))
    plt.errorbar(ra, dec, yerr=dec_err, fmt=".r", zorder=0, label="data")
    plt.plot(track_radec[:, 0], track_radec[:, 1], label="track")
    plt.scatter(ra, track_interp(ra), s=10, c="k", label="match")
    frame1.set_xticklabels([])
    frame1.set(xlim=xlim, ylim=ylim, ylabel="Dec [degree]")
    plt.grid(True)

    # residual plot
    frame2 = fig.add_axes((0.1, 0.1, 0.8, 0.2))
    plt.axhline(0, c="tab:blue")
    plt.plot(ra, residual, "or")
    frame2.set(xlim=xlim, ylim=ylim, xlabel="RA [degree]")
    plt.grid(False)

    plt.set(**kw, ax=frame1)

    return fig, (frame1, frame2)


# /def


##############################################################################
# END
