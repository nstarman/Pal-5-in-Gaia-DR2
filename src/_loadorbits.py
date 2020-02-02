# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : loadorbits
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""load orbits.

Routine Listings
----------------
load_orbit_and_stream
plot_stream_orbit_mollweide
plot_spatial_closeup
plot_spatial
plot_velocity
plot_angular_momentum
load_Pal5_orbit_and_stream

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "load_orbit_and_stream",
    "plot_stream_orbit_mollweide",
    "plot_spatial_closeup",
    "plot_spatial",
    "plot_velocity",
    "plot_angular_momentum",
    "load_Pal5_orbit_and_stream",
]


#############################################################################
# IMPORTS

# GENERAL
import numpy as np  # numerical python

# galpy
from galpy.potential import MWPotential2014
from galpy.orbit import Orbit
from galpy.util import bovy_plot, bovy_coords

# astropy
from astropy import units as u, coordinates as coord
from astropy.table import Table, QTable
from astropy.coordinates.representation import CartesianDifferential

# CUSTOM
from astroPHD.plot import starkplot as plt
from astroPHD import ObjDict, LogFile

# PROJECT-SPECIFIC
from .progenitors import loadProgenitor
from .util.coordinates.frameattrs import convert_repr
from .orbit import SequentialOrbits


##############################################################################
# SETUP

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


##############################################################################
# CODE
##############################################################################

##############################################################################
# LOAD ORBITS


def load_orbit_and_stream(
    nbodypath,
    adj=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    nbodyformat="ascii.commented_header",
    orbitpath=None,
    orbitname=None,
    plot=False,
    returnplots=False,
    logger=_LOGFILE,
    verbose=None,
):
    """Load Orbit and N-Body data.

    N-Body requires Galactocentric coordinates

    Parameters
    ----------
    nbodypath: str
        relative path to nbody data
    adj: list
        spatial and kinematic adjustments to the N-Body data
        form: [x (kpc), y (kpc), z (kpc), v_x (km/s), v_y (km/s), v_z (km/s)]
        item elements must have units compatible with listed in 'form'
            or no units, in which case they are assumed to have above unit
    orbitpath: str (or None)
        relative path to orbit .json file
        see .json form below in 'Examples'
    plot: bool
        whether to make and show plots
    returnplots: bool
        whether to return plots

    Returns
    -------
    orbit: ObjDict
        orbit properties
        {orb / orbb: galpy Orbit
            integrated forward/backward orbit
         sc / scb: astropy SkyCoord
            forward / backward SkyCoord from orbit
         gc / gcb: astropy SkyCoord
            Galactocentric coordinates
         gal / galb: astropy SkyCoord
            Galactic coordinates
         eq / eqb: astropy SkyCoord
            Equatorical coordinates
         Lgc / Lgcb: ndarray
            forward / backward angular moment array
            Galactocentric coordinates
         Lgal / Lgalb: ndarray
             forward / backward angular moment array
             Galactic coordinates
        }

    stream: ObjDict
        stream properties
        {orb: astropy QTable
            N-Body Table
         sc: astropy SkyCoord
            SkyCoord from orbit
         gc: astropy SkyCoord
            Galactocentric coordinates
         gal: astropy SkyCoord
            Galactic coordinates
         eq: astropy SkyCoord
            Equatorial coordinates
         Lgc: ndarray
            angular moment array
            Galactocentric coordinates
         Lgal: ndarray
            angular moment array
            Galactic coordinates
        }

    plot: matplotlib.pyplot figure

    """
    def sphericalcoslat(sc, wrap_angle=None):
        # sc.representation_type = 'spherical'
        # sc.differential_type = coord.SphericalCosLatDifferential
        # if wrap_angle is not None:
        #     getattr(sc, wrap_angle).wrap_angle = 180 * u.deg  # TODO

        return convert_repr(
            sc, representation="spherical", differential="sphericalcoslat"
        )

    # ------------------------------------------------------------------------
    # Progenitor

    if orbitpath is None and orbitname is None:
        raise ValueError("orbitpath and orbitname is not None")
    elif orbitpath is None:
        prog = loadProgenitor(orbitname, logger=logger)
    else:
        raise NotImplementedError("not yet implemented")

    logger.report(f"loading Progenitor {prog.name} orbit", verbose=verbose)

    # Parameter Values
    vxvv = prog.coord
    ro = 8.0 * u.kpc
    vo = 220 * u.km / u.s
    zo = 0.025 * u.kpc
    vsun = [-11.1, 24.0, 7.25] * u.km / u.s  # -11.1 b/c handedness
    vsunGC = [11.1, 244, 7.25] * u.km / u.s  # +11.1 b/c handedness

    # Making Orbit (Galactic)
    t_arr = np.linspace(0, prog["t+"]["end"], num=prog["t+"]["num"])
    t_arr = t_arr * u.Gyr if not hasattr(t_arr, "unit") else t_arr
    tb_arr = np.linspace(0, prog["t-"]["end"], num=prog["t-"]["num"])
    tb_arr = tb_arr * u.Gyr if not hasattr(tb_arr, "unit") else tb_arr

    o = SequentialOrbits(
        vxvv=vxvv, ro=ro, vo=vo, zo=zo, solarmotion=vsun, pot=MWPotential2014
    )
    o.integrate(t_arr)

    # import pdb; pdb.set_trace()

    o.add_backward_orbit()
    o.integrate(tb_arr)

    logger.report("\tmade orbit", verbose=verbose, start_at=2)

    # ---------------------------------------------------
    # SkyCoordinates

    orbSkyCrd = sphericalcoslat(o.SkyCoord("full"))
    # Galactocentric (x, y, z, v_xyz)
    gc_orb = sphericalcoslat(orbSkyCrd.galactocentric, "lon")
    # Galactic (l, b)
    gal_orb = sphericalcoslat(orbSkyCrd.galactic, "l")
    # Equatorial (ra, dec)
    eq_orb = sphericalcoslat(orbSkyCrd.icrs, "ra")

    logger.report("\tmade skycoords", verbose=verbose, start_at=2)

    # ---------------------------------------------------
    # Angular Momentum
    # Angular Momentum in the Galactocentric Reference Frame
    Lgc_orb = o.L("full")

    # # Angular Momentum in the Galactic Reference Frame
    Lgal_orb = np.cross(
        gal_orb.cartesian.xyz.transpose(), gal_orb.velocity.d_xyz.transpose()
    )

    logger.report("\tmade angular momentum", verbose=verbose, start_at=2)

    # ---------------------------------------------------
    # Storing

    orbit = ObjDict(
        "orbit",
        orb=o,
        sc=orbSkyCrd,
        # skycoords
        eq=eq_orb,
        icrs=eq_orb,
        gc=gc_orb,
        galactocentric=gc_orb,
        gal=gal_orb,
        galactic=gal_orb,
        # angular momentum
        Lgc=Lgc_orb,
        Lgal=Lgal_orb,
        # info
        info=prog,
    )

    logger.report("\tstored orbit", verbose=verbose, start_at=2)

    # ------------------------------------------------------------------------
    # N-body Stream Data

    logger.report("making N-body", verbose=verbose)

    # print(nbodypath, nbodyformat)
    df = Table.read(nbodypath, format=nbodyformat)

    # units
    logger.report("\tfixing units", verbose=verbose, start_at=2)
    udict = {
        "mass": u.Msun,
        "x": u.kpc,
        "y": u.kpc,
        "z": u.kpc,
        "v_x": u.km / u.s,
        "v_y": u.km / u.s,
        "v_z": u.km / u.s,
    }
    for key, unit in udict.items():
        if key in df.colnames:
            setattr(df[key], "unit", unit)

    # Manual Adjustments
    logger.report("\tadjusted position", verbose=verbose, start_at=2)
    for i, unit in enumerate([*[u.kpc] * 3, *[u.km / u.s] * 3]):
        try:
            adj[i].to(unit)
        except AttributeError:  # doesn't have units
            adj[i] *= unit  # assigning correct units
        except u.UnitConversionError:  # not in correct units
            raise u.UnitConversionError(
                f"{adj[i]} is not compatible with {unit}"
            )

    df["x"] = (df["x"] + adj[0]).to(u.kpc)
    df["y"] = (df["y"] + adj[1]).to(u.kpc)
    df["z"] = (df["z"] + adj[2]).to(u.kpc)
    df["v_x"] = (df["v_x"] + adj[3]).to(u.km / u.s)
    df["v_y"] = (df["v_y"] + adj[4]).to(u.km / u.s)
    df["v_z"] = (df["v_z"] + adj[5]).to(u.km / u.s)

    # ---------------------------------------------------
    # SkyCoords

    streamSkyCrd = coord.SkyCoord(
        x=-df["x"],
        y=df["y"],
        z=df["z"],
        v_x=-df["v_x"],  # -v_x b/c handedness
        v_y=df["v_y"],
        v_z=df["v_z"],
        galcen_distance=np.sqrt(ro ** 2 + zo ** 2),
        # galcen_distance=8.3 * u.kpc,
        z_sun=zo,
        galcen_v_sun=CartesianDifferential(*vsunGC),
        frame="galactocentric",
        # galcen_coord=coord.ICRS(ra=266.4051 * u.deg, dec=-28.936175 * u.deg)
    ).icrs

    streamSkyCrd = sphericalcoslat(streamSkyCrd, "ra")
    # streamSkyCrd.representation_type = 'spherical'
    # streamSkyCrd.differential_type = coord.SphericalCosLatDifferential

    # Galactocentric (x, y, z, v_xyz)
    gc_stream = sphericalcoslat(streamSkyCrd.galactocentric, "lon")

    # Galactic (l, b)
    gal_stream = sphericalcoslat(streamSkyCrd.galactic, "l")

    # Equatorial (ra, dec)
    eq_stream = sphericalcoslat(streamSkyCrd.icrs, "ra")

    logger.report("\tmade skycoords", verbose=verbose, start_at=2)

    # ---------------------------------------------------
    # Galactocentric Angular Momentum
    Lgc_stream = np.cross(
        df["x", "y", "z"].to_pandas().values,
        df["v_x", "v_y", "v_z"].to_pandas().values,
    )

    # Lgc_stream = np.cross(   # DIDN'T TRANSFORM CORRECTLY
    #     gc_stream.cartesian.xyz.transpose(),
    #     gc_stream.velocity.d_xyz.transpose())

    # Angular Momentum in the Galactic Reference Frame
    Lgal_stream = np.cross(
        gal_stream.cartesian.xyz.transpose(),
        gal_stream.velocity.d_xyz.transpose(),
    )
    # np.cross(uvw, vuvvvw)

    logger.report("\tmade angular momentum", verbose=verbose, start_at=2)

    # ---------------------------------------------------
    df = QTable(df)

    stream = ObjDict(
        "stream",
        orb=df,
        sc=streamSkyCrd,
        # skycoords
        gc=gc_stream,
        galactocentric=gc_stream,
        gal=gal_stream,
        galactic=gal_stream,
        eq=eq_stream,
        icrs=eq_stream,
        # angular momentum
        Lgc=Lgc_stream,
        Lgal=Lgal_stream,
    )

    logger.report("\tstored stream", verbose=verbose, start_at=2)

    # ------------------------------------------------------------------------
    # Plotting

    if plot is True:

        logger.report("plotting", verbose=verbose)

        fig0 = plot_stream_orbit_mollweide(eq_stream, eq_orb)
        fig1 = plot_spatial_closeup(gc_orb, gc_stream)  # closeup Plots
        fig2 = plot_spatial(gal_stream, gal_orb)  # Position
        fig3 = plot_velocity(gc_orb, gc_stream)  # Velocity
        fig4 = plot_angular_momentum(Lgal_orb, Lgal_stream)  # Angular Momentum

        plt.show()

        if returnplots is True:
            print("done loading orbits")
            return orbit, stream, (fig0, fig1, fig2, fig3, fig4)

    else:
        logger.report("not plotting", verbose=verbose)

    # ------------------------------------------------------------------------

    logger.report("done loading orbits", verbose=verbose)

    if returnplots is True:
        return orbit, stream, None
    return orbit, stream


# /def


def _scatterhelp(ax, source, q1="ra", q2="dec", label="", wrap=False, **kw):
    if wrap:
        res = ax.scatter(
            getattr(source, q1).wrap_at(180 * u.deg).rad,
            getattr(source, q2).wrap_at(180 * u.deg).rad,
            label=label,
            **kw,
        )
    else:
        res = ax.scatter(
            getattr(source, q1), getattr(source, q2), label=label, **kw
        )
    return res


# /def


# ----------------------------------------------------------------------------


def _scatterhelpradec(ax, source, label, **kw):
    return _scatterhelp(ax, source, q1="ra", q2="dec", label=label, wrap=True)


# /def


# ----------------------------------------------------------------------------


@plt.mpl_decorator(
    fig="new",
    figsize=(8, 6),
    title=r"$\alpha$ & $\delta$ (Equatorial)",
    xlabel=r"$\alpha$",
    ylabel=r"$\delta$",
)
def plot_stream_orbit_mollweide(eq_orb, eq_stream, **kw):
    """plot_stream_orbit_mollweide."""
    fig = plt.gcf()
    ax = fig.add_subplot(111, projection="mollweide")

    # stream
    l0 = _scatterhelpradec(ax, eq_stream, "stream", s=2)
    l1 = _scatterhelpradec(ax, eq_orb, "progenitor", s=2, c="r")
    l2 = _scatterhelpradec(ax, eq_orb[0], "Pal5", s=30, c="r", edgecolor="k")

    return fig


# /def


# ----------------------------------------------------------------------------


def plot_spatial_closeup(gc_orb, gc_stream):
    """plot_spatial_closeup.

    TODO
    ----
    make general

    """
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    qpairs = [("x", "y"), ("x", "z"), ("y", "z")]
    q1lims = [(8, 8.4), (8, 8.4), (0, 0.5)]
    q2lims = [(-0.5, 1), (16.2, 16.8), (16.2, 16.8)]

    for ax, (q1, q2), q1lim, q2lim in zip(axs, qpairs, q1lims, q2lims):
        # Stream
        dd1 = getattr(gc_stream.cartesian, q1).value
        dd2 = getattr(gc_stream.cartesian, q2).value
        _ind = (
            (q1lim[0] < dd1)
            & (dd1 < q1lim[1])
            & (q2lim[0] < dd2)
            & (dd2 < q2lim[1])
        )

        plt.sca(ax)
        # TODO switch to _scatterhelp
        plt.scatter(
            getattr(gc_stream.cartesian, q1)[_ind],
            getattr(gc_stream.cartesian, q2)[_ind],
            s=2,
            label="stream",
        )
        # Pal5 Integrated
        plt.scatter(
            getattr(gc_orb.cartesian, q1),
            getattr(gc_orb.cartesian, q2),
            s=1,
            c="r",
            label="Pal5",
        )
        # Progenitor
        plt.scatter(
            getattr(gc_orb.cartesian, q1)[0],
            getattr(gc_orb.cartesian, q2)[0],
            s=30,
            c="r",
            edgecolor="k",
            label="Pal5",
        )

        plt.set(
            title="{q1} & {q2} (Galactocentric)".format(q1=q1, q2=q2),
            xlabel=q1,
            ylabel=q2,
            xlim=q1lim,
            ylim=q2lim,
        )
    fig.tight_layout()
    return fig


# /def


# ----------------------------------------------------------------------------


def plot_spatial(gal_stream, gal_orb):
    """plot_spatial."""
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (q1, q2) in zip(axs, (("x", "y"), ("x", "z"), ("y", "z"))):
        s1_arr = getattr(gal_stream.cartesian, q1)[::50]
        s2_arr = getattr(gal_stream.cartesian, q2)[::50]
        ax.scatter(s1_arr, s2_arr, s=2, label="stream")

        p1_arr = getattr(gal_orb.cartesian, q1)
        p2_arr = getattr(gal_orb.cartesian, q2)
        ax.scatter(p1_arr, p2_arr, s=2, c="r", label="progenitor")

        plt.set(
            ax=ax,
            title="{q1} & {q2} (Galactic)".format(q1=q1, q2=q2),
            xlabel=r"{q1} [{unit}]".format(q1=q1, unit=s1_arr.unit),
            ylabel=r"{q2} [{unit}]".format(q2=q2, unit=s2_arr.unit),
        )
    fig.tight_layout()

    return fig


# /def


# ----------------------------------------------------------------------------


def plot_velocity(gc_orb, gc_stream):
    """Plot velocity.

    Parameters
    ----------
    gc_orb: SkyCoord
        orbit in galactocentric coordinates
    gc_stream: SkyCoord
        stream in galactocentric coordinates

    Return
    ------
    fig: pyplot.Figure

    """
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    vpairs = (("d_x", "d_y"), ("d_x", "d_z"), ("d_y", "d_z"))
    v1lims = [(-100, -20), (-100, -20), (-130, -100)]
    v2lims = [(-130, -100), (-100, 25), (-100, 25)]
    streamvals = gc_stream.cartesian.differentials["s"]
    orbvals = gc_orb.cartesian.differentials["s"]

    for ax, (q1, q2), q1lim, q2lim in zip(axs, vpairs, v1lims, v2lims):

        dd1 = getattr(streamvals, q1).value
        dd2 = getattr(streamvals, q2).value
        _ind = (
            (q1lim[0] < dd1)
            & (dd1 < q1lim[1])
            & (q2lim[0] < dd2)
            & (dd2 < q2lim[1])
        )

        s1_arr = getattr(streamvals, q1)[_ind]
        s2_arr = getattr(streamvals, q2)[_ind]
        ax.scatter(s1_arr, s2_arr, s=2, label="stream")

        p1_arr = getattr(orbvals, q1)
        p2_arr = getattr(orbvals, q2)
        ax.scatter(p1_arr, p2_arr, s=2, c="r", label="progenitor")

        # progenitor
        ax.scatter(
            getattr(orbvals, q1)[0],
            getattr(orbvals, q2)[0],
            s=30,
            c="r",
            edgecolor="k",
            label="Pal5",
        )

        # plot properties
        ax.set_title("{q1} & {q2} (Galactocentric)".format(q1=q1, q2=q2))
        ax.set_xlabel(r"{q1} [{unit}]".format(q1=q1, unit=s1_arr.unit))
        ax.set_ylabel(r"{q2} [{unit}]".format(q2=q2, unit=s2_arr.unit))
        ax.set_xlim(q1lim)
        ax.set_ylim(q2lim)
        ax.legend()

    # /for

    fig.tight_layout()

    return fig


# /def


# ----------------------------------------------------------------------------


def plot_angular_momentum(Lgal_orb, Lgal_stream):
    """plot_angular_momentum.

    Parameters
    ----------
    Lgal_orb: ndarray
    Lgal_stream: ndarray

    Returns
    -------
    fig: pyplot.Figure

    """
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    for ax, (q1, q2) in zip(axs, (("x", "y"), ("x", "z"), ("y", "z"))):

        opts = ["x", "y", "z"]
        i, j = opts.index(q1), opts.index(q2)

        ax.scatter(Lgal_stream[:, i], Lgal_stream[:, j], s=2, label="stream")

        ax.scatter(
            Lgal_orb[:, i], Lgal_orb[:, j], s=2, c="r", label="progenitor"
        )

        ax.set_title("$L_{q1}$ & $L_{q2}$ (Galactic)".format(q1=q1, q2=q2))
        ax.set_xlabel(r"$L_{q1}$ [{unit}]".format(q1=q1, unit="todo"))
        ax.set_ylabel(r"$L_{q2}$ [{unit}]".format(q2=q2, unit="todo"))
        ax.legend()

    # /for

    fig.tight_layout()

    return fig


# /def


##############################################################################
# loadPal5orbits


def load_Pal5_orbit_and_stream(
    nbodypath,
    adj=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    nbodyformat="ascii.ecsv",
    plot=False,
    returnplots=False,
    logger=_LOGFILE,
    verbose=None,
):
    """Load Pal 5 Orbit and N-Body data.

    Parameters
    ----------
    nbodypath: str
        relative path to nbody data
    adj: list
        spatial and kinematic adjustments to the N-Body data
        form: [x (kpc), y (kpc), z (kpc), v_x (km/s), v_y (km/s), v_z (km/s)]
        item elements must have units compatible with listed in 'form'
            or no units, in which case they are assumed to have above unit
    plot: bool
        whether to make and show plots
    returnplots: bool
        whether to return plots

    Returns
    -------
    pal5 : dict
    stream : dict
    plot : matplotlib.pyplot.Figure

    """
    return load_orbit_and_stream(
        nbodypath,
        adj=adj,
        nbodyformat=nbodyformat,
        orbitname="Palomar 5",
        plot=plot,
        returnplots=returnplots,
        logger=logger,
        verbose=verbose,
    )


# /def

##############################################################################
# End
