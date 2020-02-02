# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : MegaCamGen1_from_PS1
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""MegaCamGen1_from_PS1."""

__author__ = "Nathaniel Starkman"
__credits__ = [
    "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html"
]


##############################################################################
# IMPORTS

# GENERAL
import warnings

from astropy import units
from astropy.table import Table

# PROJECT-SPECIFIC
from .units_decorators import quantity_io


#############################################################################
# CODE
#############################################################################


@quantity_io()
def U_MP9301(ps: Table, **kw) -> units.mag:
    """U_MP9301.

    gmi = (gPS-iPS)
    uCFHT - gPS = .523 - .343 gmi + 2.44 gmi^2 - .998 gmi^3

    limited to .3 mag < gmi < 1.5 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. Top row, 1st plot

    Parameters
    ----------
    ps: Astropy Table
        need: g col
        either: i, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    U_MP9301: array_like
        CFHT u-band

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (0.3 * units.mag < gmi) & (gmi < 1.5 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.U: not all .3 mag < (g-i)_ps < 1.5 mag")

    c0 = 0.523 * units.mag
    c1 = -0.343
    c2 = 2.44 / units.mag
    c3 = -0.998 / units.mag ** 2
    g_ps = ps[g]

    u_cfht = g_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return u_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def G_MP9401(ps: Table, **kw) -> units.mag:
    """G_MP9401.

    gmi = (gPS-iPS)
    gCFHT - gPS = -.001 - .004 gmi - .0056 gmi^2 + .00292 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. Top row, 2nd plot

    Parameters
    ----------
    ps: Astropy Table
        need: g col
        either: i, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    G_MP9401: array_like
        CFHT g-band

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.G: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = -0.001 * units.mag
    c1 = -0.004
    c2 = -0.0056 / units.mag
    c3 = 0.00292 / units.mag ** 2
    g_ps = ps[g]

    g_cfht = g_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return g_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def R_MP9601(ps: Table, **kw) -> units.mag:
    """R_MP9601.

    gmi = (gPS-iPS)
    rCFHT - rPS = .002 - .017 gmi + .00554 gmi^2 - .000692 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. Top row, 3rd plot

    Parameters
    ----------
    ps: Astropy Table
        need: r, col
        either: (g & i), g-i col

    Other Parameters
    ----------------
    g: str   (default 'g')
        g column name
    r: str   (default 'r')
        r column name
    i: str   (default 'i')
        i column name
    gmi: str   (default 'g-i')
        g-i column name

    Returns
    -------
    r_cfht: array_like

    """
    g, r, i = kw.get("g", "g"), kw.get("r", "r"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.R: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = 0.002 * units.mag
    c1 = -0.017
    c2 = 0.00554 / units.mag
    c3 = -0.000692 / units.mag ** 2
    r_ps = ps[r]

    r_cfht = r_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return r_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def I_MP9701(ps: Table, **kw) -> units.mag:
    """I_MP9701.

    gmi = (gPS-iPS)
    iCFHT - iPS = .001 - .021 gmi + .00398 gmi^2 - .00369 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 2nd row, 1st plot

    Parameters
    ----------
    ps: Astropy Table
        need: i col
        either: g, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    i_cfht: array_like

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.I: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = 0.001 * units.mag
    c1 = -0.021
    c2 = 0.00398 / units.mag
    c3 = -0.00369 / units.mag ** 2
    i_ps = ps[i]

    i_cfht = i_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return i_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def Z_MP9801(ps: Table, **kw) -> units.mag:
    """Z_MP9801.

    gmi = (gPS-iPS)
    zCFHT - zPS = -.009 - .029 gmi + .012 gmi^2 - .00367 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 2nd row, 2nd plot

    Parameters
    ----------
    ps: Astropy Table
        need: z col
        either: (g & i), g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    z: str  (default 'z')
        z column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    z_cfht: array_like

    """
    g, i, z = kw.get("g", "g"), kw.get("i", "i"), kw.get("z", "z")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.Z: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = -0.009 * units.mag
    c1 = -0.029
    c2 = 0.012 / units.mag
    c3 = -0.00367 / units.mag ** 2
    z_ps = ps[z]

    z_cfht = z_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return z_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def I_MP9702(ps: Table, **kw) -> units.mag:
    """I_MP9702.

    gmi = (gPS-iPS)
    iCFHT - iPS = -.005 + .004 gmi + .0124 gmi^2 - .0048 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 2nd row, 3rd plot

    Parameters
    ----------
    ps: Astropy Table
        need: i col
        either: g, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    z_cfht: array_like

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.I: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = -0.005 * units.mag
    c1 = +0.004
    c2 = 0.0124 / units.mag
    c3 = -0.0048 / units.mag ** 2
    i_ps = ps[i]

    z_cfht = i_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return z_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def U_MP9302(ps: Table, **kw) -> units.mag:
    """I_MP9702.

    gmi = (gPS-iPS)
    uCFHT - gPS = .823 - 1.36 gmi + 4.18 gmi^2 - 1.64 gmi^3

    limited to .3 mag < gmi < 1.5 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 3rd row, 1st plot

    Parameters
    ----------
    ps: Astropy Table
        need: g col
        either: i, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    u_cfht: array_like

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (0.3 * units.mag < gmi) & (gmi < 1.5 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.U: not all .3 mag < (g-i)_ps < 1.5 mag")

    c0 = 0.823 * units.mag
    c1 = -1.360
    c2 = 4.18 / units.mag
    c3 = -1.64 / units.mag ** 2
    g_ps = ps[g]

    u_cfht = g_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return u_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def G_MP9402(ps: Table, **kw) -> units.mag:
    """G_MP9402.

    gmi = (gPS-iPS)
    gCFHT - gPS = .014 - .059 gmi - .00313 gmi^2 - .00178 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 3rd row, 2nd plot

    Parameters
    ----------
    ps: Astropy Table
        need: g col
        either: i, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    g_cfht: array_like

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.G: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = 0.014 * units.mag
    c1 = 0.059
    c2 = -0.00313 / units.mag
    c3 = -0.00178 / units.mag ** 2
    g_ps = ps[g]

    g_cfht = g_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return g_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def R_MP9602(ps: Table, **kw) -> units.mag:
    """R_MP9602.

    gmi = (gPS-iPS)
    rCFHT - rPS = .003 - .05 gmi - .0125 gmi^2 - .00699 gmi^3

    limited to -1 mag < gmi < 3 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 3rd row, 3rd plot

    Parameters
    ----------
    ps: Astropy Table
        need: r col
        either: (g & i), g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    r: str  (default 'r')
        r column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    r_cfht: array_like

    """
    g, r, i = kw.get("g", "g"), kw.get("r", "r"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 3.0 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.R: not all -1 mag < (g-i)_ps < 3 mag")

    c0 = 0.003 * units.mag
    c1 = -0.050
    c2 = 0.0125 / units.mag
    c3 = -0.00699 / units.mag ** 2
    r_ps = ps[r]

    r_cfht = r_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return r_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def I_MP9703(ps: Table, **kw) -> units.mag:
    """I_MP9703.

    gmi = (gPS-iPS)
    iCFHT - iPS = .006 - .024 gmi + .00627 gmi^2 - .00523 gmi^3

    limited to -1 mag < gmi < 3.6 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 4th row, 1st plot

    Parameters
    ----------
    ps: Astropy Table
        need: i col
        either: g, g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    i_cfht: array_like

    """
    g, i = kw.get("g", "g"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 3.6 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.I: not all -1 mag < (g-i)_ps < 3.6 mag")

    c0 = 0.006 * units.mag
    c1 = -0.024
    c2 = 0.00627 / units.mag
    c3 = -0.00523 / units.mag ** 2
    i_ps = ps[i]

    i_cfht = i_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return i_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def Z_MP9901(ps: Table, **kw) -> units.mag:
    """Z_MP9901.

    gmi = (gPS-iPS)
    zCFHT - zPS = -.016 - .069 gmi + .0239 gmi^2 - .0056 gmi^3

    limited to -1 mag < gmi < 4 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 4th row, 2nd plot

    Parameters
    ----------
    ps: Astropy Table
        need: z col
        either: (g & i), g-i col

    Other Parameters
    ----------------
    i: str  (default 'i')
        i column name
    z: str  (default 'z')
        z column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    z_cfht: array_like

    """
    g, i, z = kw.get("g", "g"), kw.get("i", "i"), kw.get("z", "z")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 4.0 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.Z: not all -1 mag < (g-i)_ps < 4 mag")

    c0 = -0.016 * units.mag
    c1 = -0.069
    c2 = 0.0239 / units.mag
    c3 = -0.0056 / units.mag ** 2
    z_ps = ps[z]

    z_cfht = z_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return z_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def GRI_MP9605(ps: Table, **kw) -> units.mag:
    """GRI_MP9605.

    gmi = (gPS-iPS)
    gCFHT - rPS = .005 + .244 gmi - .0692 gmi^2 - .0014 gmi^3

    limited to -1 mag < gmi < 1.2 mag

    filter transformations from
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    Pan-STARRS to MegaCam plots. 4th row, 3rd plot

    Parameters
    ----------
    ps: Astropy Table
        need: r col
        either: (g & i), g-i col

    Other Parameters
    ----------------
    g: str  (default 'g')
        g column name
    r: str  (default 'r')
        r column name
    i: str  (default 'i')
        i column name
    gmi: str  (default 'g-i')
        g-i column name

    Returns
    -------
    r_cfht: array_like

    """
    g, r, i = kw.get("g", "g"), kw.get("r", "r"), kw.get("i", "i")
    gmi = kw.get("gmi", "g-i")

    if gmi in ps.colnames:
        gmi = ps[gmi]
    else:
        gmi = ps[g] - ps[i]

    ind = (-1.0 * units.mag < gmi) & (gmi < 1.2 * units.mag)
    if not all(ind):
        warnings.warn("MCg1.G: not all -1 mag < (g-i)_ps < 1.2 mag")

    c0 = 0.005 * units.mag
    c1 = 0.244
    c2 = -0.0692 / units.mag
    c3 = -0.0014 / units.mag ** 2
    r_ps = ps[r]

    r_cfht = r_ps + c0 + (c1 * gmi) + (c2 * gmi ** 2) + (c3 * gmi ** 3)
    return r_cfht


# /def


# -------------------------------------------------------------------------


@quantity_io()
def GmR(ps: Table, gfilter="9401", rfilter="9601", **kw) -> units.mag:
    """GmR.

    Parameters
    ----------
    gfilter: {"9401", "9402"}
    rfilter: {"9601", "9602"}

    Returns
    -------
    G: array_like

    """
    if gfilter in ("9401", "9402"):
        G = globals().get("G_MP" + gfilter)
    else:
        gopts = ("9401", "9402")
        raise ValueError(f"{gfilter} wrong not in {gopts}")

    if rfilter in ("9601", "9602"):
        R = globals().get("R_MP" + rfilter)
    else:
        ropts = ("9601", "9602")
        raise ValueError(f"{rfilter} wrong not in {ropts}")

    return G(ps, **kw) - R(ps, **kw)


# /def
