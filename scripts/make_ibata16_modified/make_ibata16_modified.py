#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : make Ibata 2016 CFHT table
# AUTHOR  : Nathaniel Starkman
# PROJECT : Palomar 5 in Gaia DR2
#
# ----------------------------------------------------------------------------

### Docstring and Metadata
"""script: Make Modified Ibata 2016 Table


"""

__author__ = "Nathaniel Starkman"

##############################################################################
### Imports

## General
import os, sys, time, pdb, copy
# import numpy as np

## Astropy
from astropy.table import Table, QTable

## Gaia-tools
from gaia_tools.util import table_utils as tutil

## Custom
sys.path.insert(0, '../../')
from src.util.logging import LogFile

## Project-Specific


##############################################################################
### Parameters & Setup

# General
_PLOT = True                                # Plot the output

# Data
_IN_FILE = '../../data/ibata2016cfht/pal5.fit'
_IN_FILE_FORMAT = 'fits'

_OUT_FILE = '../../data/ibata2016cfht/pal5_modified.fit'
_OUT_FILE_FORMAT = 'fits'

# Log file
print(__name__)
_LOGFILE = LogFile.open('./log.txt', verbose=2, mode='w',
                        header='Make Ibata 2016 Modified')

# ----------------------------------------------------------------------------
### Setup

_LOGFILE.write('opening ' + _IN_FILE)
dfo = Table.read(_IN_FILE, format=_IN_FILE_FORMAT);

_LOGFILE.write(f'opened\ntable columns: {dfo.colnames}', endsection=True)

##############################################################################
### Running the Script

# Renaming columns
_LOGFILE.write('renaming columns to...')
tutil.rename_columns(
    dfo,
    ['E_B-V_', 'E(B-V)'],
    _RAJ2000='ra', _DEJ2000='dec',
    g0mag='g', e_g0mag='g_err', r0mag="r", e_r0mag='r_err',
    _Glon='L', _Glat='B'
)
_LOGFILE.write(dfo.colnames)

# Adding color column
_LOGFILE.write('adding g-r color column', start='\n')
tutil.add_color_col(dfo, 'g', 'r', color='g-r')

# New Table
_LOGFILE.write('Making new QTable w/ all columns except RAJ2000, DEJ2000',
               start='\n')
names = tutil.drop_colnames(dfo.colnames, 'RAJ2000', 'DEJ2000')
df = QTable(dfo[names])
_LOGFILE.write(f'new table columns: {df.colnames}', start='\n')

# Writing
_LOGFILE.write('saving to ' + _OUT_FILE, start='\n')
df.write(_OUT_FILE, format=_OUT_FILE_FORMAT, overwrite=True)

##############################################################################
### Closing
_LOGFILE.close()
