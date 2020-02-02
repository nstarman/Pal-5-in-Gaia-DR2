# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : src initialization
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""AST1501.

Routine Listings
----------------
makefilepaths
load_orbit_and_stream
load_Pal5_orbit_and_stream

"""

__author__ = "Nathaniel Starkman"

##############################################################################
# IMPORTS

# CUSTOM
from astroPHD import LogFile, ObjDict

# PROJECT-SPECIFIC
from ._loadorbits import load_orbit_and_stream, load_Pal5_orbit_and_stream


###############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


###############################################################################
# CODE


def makefilepaths(
    datapath="/geir_data/scr/nstarkman/Pal-5-in-Gaia-DR2/",
    nbody="data/nbody/pal5.dat",
    nbodyformat="ascii.ecsv",
    skydata="scripts/get_gaia_window/base/output/window.fits",
    logger=_LOGFILE,
    verbose=None,
    **kw,
):
    """Make Notebook Options.

    INPUT
    -----
    datapath: str  (default = '/geir_data/scr/nstarkman/data/')
        the file path to the data
    nbody: str  (default = 'nbody3')
        what nbody to use
        f'{datapath}{nbody}/nbody/pal5.dat'
    nbodyformat: str  (default = 'ascii.ecsv')
        what format the nbody is
    skydata: str  (default =  '040221')
        the folder in which the sky window is located
        f'{datapath}{nbody}/skydata/{skydata}_{modified}.fits'
    modified: str  (default = 'Modified')
        modifier for skydata
        options: None, 'Modified', 'Modified_pmsub'

    OUTPUT
    ------
    opts: ObjDict
        return order:
            'nbody', 'skydata',
            'datapath', 'nbodypath', 'skydatapath',
            'nbodyformat', 'use_modified',

    """
    # the path of the N-Body
    nbodypath = datapath + nbody

    # the data path
    # if modified in ('Modified', 'Modified_pmsub'):
    #     skydatapath = f'{datapath}/{nbody}/skydata/{skydata}_{modified}.fits'
    # else:
    #     skydatapath = f'{datapath}/{nbody}/skydata/{skydata}.fits'
    skydatapath = datapath + skydata

    opts = ObjDict(
        "options",
        datapath=datapath,
        nbody=nbody,
        nbodyformat=nbodyformat,
        skydata=skydata,
        nbodypath=nbodypath,
        skydatapath=skydatapath,
        **kw,
    )

    logger.write(
        f"""makefilepaths:
    loading data from {datapath}
    nbodypath:   {nbodypath}
    skydatapath: {skydatapath}
    """
    )

    return opts


# /def

# -------------------------------------------------------------------------
