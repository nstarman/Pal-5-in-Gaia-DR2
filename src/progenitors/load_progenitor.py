# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : load_progenitor
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Load Progenitors."""

__author__ = "Nathaniel Starkman"


##############################################################################
# IMPORTS

# GENERAL
import os
import json
import numpy as np
import importlib

# Astropy
from astropy import units as u
from astropy.coordinates import SkyCoord

# CUSTOM
from astroPHD import LogFile, ObjDict

# PROJECT-SPECIFIC
from ..util import astrarray


##############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


##############################################################################
# Seeing Progenitors


def available_progenitors(_print=True):
    """Print list of available progenitors."""
    lookupdir = get_lookupdir()

    lookup = {}
    for k, v in lookupdir.items():
        lookup.setdefault(v, []).append(k)

    if _print:
        print(lookup)

    return lookup


# /def


# ----------------------------------------------------------------------------


def get_lookupdir():
    """Return lookup directory."""
    # lookup directory
    basedir = os.path.dirname(os.path.realpath(__file__))
    fpath = basedir + "/lookup.json"

    with open(fpath, "r") as file:
        lookup = json.load(file)

    return lookup


# /def


# ----------------------------------------------------------------------------


def available_coords_in_progenitor(name):
    """Print list of available progenitors."""
    # if no fpath is directly provided
    lookup = get_lookupdir()

    # getting actual file path from lookup
    fpath = ".{}".format(lookup[name])
    pkg = os.path.dirname(__file__)
    # delete bad possible prefixes
    pkg = pkg.replace("/", ".")  # changing to importlib format
    pkg = pkg.replace("...", "").replace("..", "")  # stripping prefixes

    # open file with information
    file = importlib.import_module(fpath, pkg)

    return file.progenitor.coords.keys()


# /def


##############################################################################
### Load Progenitors


def loadProgenitor(
    name=None, pmnsigma=7, coord="main", logger=_LOGFILE, verbose=None
):
    r"""Load progenitor from name.

    Parameters
    ----------
    name: str, None  (default None)
        name of progenitor
    pmnsigma: float
        for making proper motion ranges
    coord: str  (default = 'coord')
        key of coord to use in coord
        this is only used if 'coord' is not a SkyCoord itself
        ex:  for coord='name'
            'coord': {'name': SkyCoord(
                ra=229.018 * u.deg, dec=-0.124 * u.deg, distance=23.2 * u.kpc,
                pm_ra_cosdec=-2.296 * u.mas/u.yr, pm_dec=-2.257 * u.mas/u.yr,
                radial_velocity=-58.7 * u.km/u.s)}

    """
    logger.report(f"Loaded Progenitor {name} Info:", verbose=verbose)

    # +---------------------------------------------------+
    logger.report(
        "loaded progenitor file", verbose=verbose, start_at=2, start="\t"
    )

    # progenitor name file path directory
    lookup = get_lookupdir()

    # getting actual file path from lookup
    pkg = f"{os.path.dirname(__file__)}/{lookup[name]}.py"

    # open file with information
    spec = importlib.util.spec_from_file_location(lookup[name], pkg)
    file = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(file)
    prog = file.progenitor

    # +---------------------------------------------------+
    logger.report(
        "loaded coordinates", verbose=verbose, start_at=2, start="\t"
    )

    # coord
    if "coord" in prog:
        pass
    else:  # it's in coords
        crd = prog["coords"]
        if isinstance(crd, SkyCoord):
            prog.coord = crd
        else:
            prog.coord = crd[coord]  # kw.get('coord', 'coord')

    # parallax
    prog.prlx = prog.coord.distance.to(u.mas, equivalencies=u.parallax())

    # +---------------------------------------------------+
    # Proper Motion
    logger.report(
        "made proper motions", verbose=verbose, start_at=2, start="\t"
    )

    # RA
    errange = np.array([-1, 1]) * prog.info["pmra_err"] * pmnsigma
    prog.pm_ra_cosdec_range = prog.info["pmra"] + errange
    # Dec
    errange = np.array([-1, 1]) * prog.info["pmdec_err"] * pmnsigma
    prog.pm_dec_range = prog.info["pmdec"] + errange

    return prog


# /def


##############################################################################
# Specific Progenitors


def loadPal5(pmnsigma=7, coord="main", logger=_LOGFILE, verbose=None):
    """Load PAL 5.

    Cluster center at: `15 16 05.30 -00 06 41.0` (https://arxiv.org/abs/astro-ph/0511128)
    Approximate angular size (data flag of D in A-E):  `13.18 11.48 arcmin`  (SIMBAD)

    """
    return loadProgenitor(
        "Palomar5",
        pmnsigma=pmnsigma,
        coord=coord,
        logger=logger,
        verbose=verbose,
    )


# /def


##############################################################################
# END
