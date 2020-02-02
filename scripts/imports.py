# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : Standard Import File
# AUTHOR  : Nathaniel Starkman
# PROJECT : Pal-5-in-Gaia-DR2
#
# ----------------------------------------------------------------------------

# Docstring and MetaData
"""Standard Imports File.

Returns
-------
scipy.stats.binned_statistic, poisson
src:
    makefilepaths, load_Pal5_orbit_and_stream, progenitors, select
    .plot.plot_sky_window, plot_proper_motions, plot_data_along_orbit,
          plot_data_along_orbit_in_window, plot_data_orbit_projection
    projection.split_arms, get_data_along_orbit, get_data_along_orbit,
               select_stars_in_an_arm, digitize_star_along_stream
    select.inRange
    util.pickle, quadrature
        .isochrone.CMD3p1TableRead
        .coordinates.make_stream_frame
                    .CustomCoords.register_stream_frame

"""

__author__ = "Nathaniel Starkman"


##############################################################################
# HELPER FUNCTIONS

from astroPHD.config import __config__
from astroPHD.decorators.docstring import (
    _set_docstring_import_file_helper,
    _import_file_docstring_helper
)


##############################################################################
# IMPORTS

# GENERAL
import copy
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')

# numpy
from scipy.stats import binned_statistic, poisson


# CUSTOM
from astroPHD.plot import starkplot as plt


# PROJECT-SPECIFIC
from src import (
    makefilepaths,
    load_Pal5_orbit_and_stream,
    progenitors,
    select)

from src.orbit import SequentialOrbits

from src.select import inRange

from src.util.isochrone import CMD3p1TableRead
from src.util.coordinates import make_stream_frame
from src.util import pickle
from src.util.coordinates.CustomCoords import register_stream_frame

# from src.projection import split_arms, get_data_along_orbit
from src.projection import get_data_along_orbit
from src.projection import select_stars_in_an_arm, digitize_star_along_stream

from src.plot import (
    plot_sky_window, plot_proper_motions, plot_data_along_orbit,
    plot_data_along_orbit_in_window, plot_data_orbit_projection
)


##############################################################################
# Printing Information

@_set_docstring_import_file_helper('base', __doc__)  # doc from __doc__
def base_imports_help():
    """Help for Matplotlib base imports."""
    _import_file_docstring_helper(base_imports_help.__doc__)  # formatting
# /def


if __config__.getboolean('verbosity', 'verbose-imports'):
    base_imports_help()


##############################################################################
# END
