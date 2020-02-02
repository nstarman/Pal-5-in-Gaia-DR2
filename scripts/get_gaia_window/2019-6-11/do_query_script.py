#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : make gaia query
# AUTHOR  : Nathaniel Starkman
# PROJECT : Palomar 5 in Gaia DR2
#
# ----------------------------------------------------------------------------

### Docstring and Metadata
"""script: Make Gaia Query
"""

__author__ = "Nathaniel Starkman"

##############################################################################
### IMPORTS

# ## General
import numpy as np
# from scipy.linalg import norm
# # astropy
from astropy import units as u
# from astropy.table import Table, QTable, vstack
from astropy.coordinates import SkyCoord
# # galpy & gaia-tools & mwdust
from galpy.potential import MWPotential2014

# from gaia_tools.query import cache, make_query
# from gaia_tools.util import table_utils, json
# from gaia_tools.xmatch import xmatch

# ## Custom
# from matplotlib import pyplot as plt

import sys; sys.path.insert(0, '../../../')
from src import LogFile, ObjDict

# ## Project-Specific
from src import makefilepaths  # , load_Pal5_orbit_and_stream
from src.orbit import SequentialOrbits
from scripts.get_gaia_window.make_gaia_query import query_along_orbit
# # util
# from src.util import MegaCamGen1_from_PS1 as MCG1_PS1
# from src.util.mwdust_util import load_dust_gri
# # coordinates
# from src.util.coordinates.CustomCoords import make_stream_frame


##############################################################################
### PARAMETERS & SETUP

# TODO have this in a __name__ == '__main__':

# Log file
_VERBOSE = 2
_LOGFILE = LogFile.open(f'./log.txt', verbose=_VERBOSE,
                        header='Make Gaia Query')

## General
_PLOT = True          # Plot the output
_LOCAL = True
_DUSTSTR = 'output/ps1dust_{}.dat'

## Orbit Parameters
vxvv_lead = SkyCoord(
    ra=229.022 * u.deg, dec=-.223 * u.deg, distance=22.5 * u.kpc,
    pm_ra=-2.128 * u.mas / u.yr, pm_dec=-2.18 * u.mas / u.yr,
    radial_velocity=-50 * u.km / u.s,
    differential_type='spherical'
)

vxvv_trail = SkyCoord(
    ra=229.022*u.deg, dec=-.223*u.deg, distance=22.5*u.kpc,
    pm_ra=-2.0*u.mas/u.yr, pm_dec=-2.18*u.mas/u.yr,
    radial_velocity=-50*u.km/u.s,
    differential_type='spherical'
)

t_arr = np.linspace(0, 90, num=10000) * u.Myr

## Query Parameters
_RANDOM_INDEX = None
_PHI1MIN, _PHI1MAX = -10 * u.deg, 10 * u.deg
_PHI2MIN, _PHI2MAX = -5 * u.deg, 5 * u.deg
_PHILIMS = [_PHI1MIN, _PHI1MAX, _PHI2MIN, _PHI2MAX]

_POINTS_INDS = [0, 3000, 6000, 13500, 17500, -1]

# _PHI1MIN, _PHI1MAX = -1 * u.deg, 1 * u.deg           # just for testing
# _PHI2MIN, _PHI2MAX = -1 * u.deg, 1 * u.deg           # just for testing
# _PHILIMS = [_PHI1MIN, _PHI1MAX, _PHI2MIN, _PHI2MAX]  # just for testing
# _POINTS_INDS = [0, 3000, 6000, 13500, 17500, -1]     # just for testing

# Logging
_LOGFILE.write('Parameters:')
_LOGFILE.newsection(title='General:', div='.')
_LOGFILE.write(f'plot: {_PLOT}', 'local: {_LOCAL}', sep='\n', end='\n\n')
_LOGFILE.newsection(title='Orbit Parameters', div='.')
_LOGFILE.write('vxvv_lead:', vxvv_lead)
_LOGFILE.write('vxvv_trail:', vxvv_trail)
_LOGFILE.write(f'integrating over {t_arr[0]}:{t_arr[-1]}:{len(t_arr)}')
_LOGFILE.newsection(title='Query Parameters', div='.')
_LOGFILE.write('Phi window dimensions:', f'phi1: {_PHI1MIN} : {_PHI1MAX}',
               f'phi2: {_PHI2MIN} : {_PHI2MAX}',
               f'random_index : {_RANDOM_INDEX}',
               sep='\n\t', end='\n\n')

# ----------------------------------------------------------------------------
## Setup

_LOGFILE.newsection(title='Setup:')

opts = makefilepaths(
    datapath='../../../',
    nbody='data/nbody/pal5.dat', nbodyformat='ascii.ecsv',
    skydata='scripts/get_gaia_window/base/output/window.fits',
    duststr=_DUSTSTR,
    logger=_LOGFILE, verbose=None)

# # loading orbit
# pal5, stream = load_Pal5_orbit_and_stream(
#     nbodypath=opts.nbodypath, nbodyformat=opts.nbodyformat,
#     plot=False, returnplots=False,
#     logger=_LOGFILE
# )

_LOGFILE.newsection(title='Stream Hand Fit Orbit:', div='.')

o = SequentialOrbits(vxvv_lead).integrate(t_arr, MWPotential2014)
o.add_backward_orbit(vxvv_trail).integrate(t_arr, MWPotential2014)

stream_hand_fit = ObjDict(
    'Stream Hand Fit Orbit',
    o=o,
    icrs=o.SkyCoord(t='full', frame='icrs')
)

##############################################################################
### CODE

df, MasterFrame, crd_xfm = query_along_orbit(
    stream_hand_fit.icrs, stream_hand_fit.o[0].SkyCoord(),
    _PHILIMS,
    points_inds=_POINTS_INDS,
    _random_index=_RANDOM_INDEX,
    duststr='output/ps1dust_{}.dat',
    local=_LOCAL,
    pprint=True,
    save_window=True, save_frame=True,
    logger=_LOGFILE
)


##############################################################################
### PLOTTING

# fig, ax = plt.subplots(1, 1)
# ax.scatter(df['ra'], df['dec'])
# fig.savefig('output/radec.png')

##############################################################################
### CLOSING

stream_hand_fit.save('output/stream_hand_fit.pkl')

_LOGFILE.close()

##############################################################################
### DONE
