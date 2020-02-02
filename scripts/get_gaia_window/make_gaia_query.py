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

## General
import numpy as np
from scipy.linalg import norm
import copy
# astropy
from astropy import units as u
from astropy.table import Table, QTable, vstack
from astropy.coordinates import SkyCoord
# galpy & gaia-tools & mwdust
from galpy.potential import MWPotential2014

from gaia_tools.query import cache, make_query
from gaia_tools.util import table_utils, json
from gaia_tools.xmatch import xmatch

## Custom
import starkplot as plt

import sys; sys.path.insert(0, '../../')
from src import LogFile, ObjDict

## Project-Specific
from src import makefilepaths, load_Pal5_orbit_and_stream
from src.orbit import SequentialOrbits
# util
from src.util import MegaCamGen1_from_PS1 as MCG1_PS1
from src.util.mwdust_util import load_dust_gri
from src.util.pickle import dump
# coordinates
from src.util.coordinates.CustomCoords import make_stream_frame, register_stream_frame


##############################################################################
### PARAMETERS & SETUP

# Log file
_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible


##############################################################################
### FUNCTION

def query_at_point(orbit, point, philims,
                   dt=2 * u.Myr, is_nbody=False,
                   max_separation=.5 * u.kpc, return_transforms=False,
                   as_qtable=True,
                   # for make_query
                   local=False,
                   _random_index=None,
                   # logging
                   plot=False, pprint=False,
                   save_window=False, save_frame=False,
                   logger=_LOGFILE, verbose=None):
    """

    Parameters
    ----------
    orbit : SkyCoord
        icrs orbit

    point : SkyCoord

    philims : lists
        [phi1min, phi1max, phi2min, phi2max]

    dt : Quantity
        physical_motion timestep

    is_nbody : bool

    max_separation : Quantity

    return_transforms : bool

    local : bool

    TODO options for make_query
    """

    # -------------------------------
    # Rotated frame

    Frame, coordtransforms= make_stream_frame(
        orbit, point, method='path', dt=dt,
        is_nbody=is_nbody, max_separation=max_separation,
        return_transforms=True,
        logger=logfile, verbose=verbose)

    Frame = copy.copy(Frame)                       # TODO copying
    coordtransforms = copy.copy(coordtransforms)   # TODO copying

    rotm = coordtransforms['R_icrs_to_cust']
    ra_ngp, dec_ngp = coordtransforms['ra_ngp', 'dec_ngp']

    # -------------------------------
    # the query

    (phi1min, phi1max, phi2min, phi2max) = philims

    # user asdict
    l3userasdict = {
        'K00': rotm[0, 0], 'K01': rotm[0, 1], 'K02': rotm[0, 2],
        'K10': rotm[1, 0], 'K11': rotm[1, 1], 'K12': rotm[1, 2],
        'K20': rotm[2, 0], 'K21': rotm[2, 1], 'K22': rotm[2, 2],
        'sindecngp': np.sin(dec_ngp),
        'cosdecngp': np.cos(dec_ngp),
        'mcosdecngp': -np.cos(dec_ngp),
        'mrangp': -ra_ngp.to_value('deg'),
        'phi1min': phi1min.to_value('rad'), 'phi1max': phi1max.to_value('rad'),
        'phi2min': phi2min.to_value('rad'), 'phi2max': phi2max.to_value('rad'),
        # 'pm_phi1min': pmphi1min, 'pm_phi1max': pmphi1max,
        # 'pm_phi2min': pmphi2min, 'pm_phi2max': pmphi2max
    }

    with open('skydataselection.json', 'r') as file:
            querystr = json.strjoinall(json.load(file))

    df = make_query(
        # gaia_mags=True,
        panstarrs1=True, defaults='full',
        user_cols=querystr['l3cols'], use_AS=True, user_ASdict=l3userasdict,
        WHERE=querystr['l3sel'],
        ORDERBY='gaia.source_id',
        pprint=pprint,
        # doing query
        do_query=True, units=True, local=local,
        # cache=f'pal5_{opts.nbody}_phi1_m25p25_phi2_m5p20',  # TODO

        # Inner Query
        FROM=make_query(
            # gaia_mags=True,
            user_cols=querystr['l2cols'], defaults='full',
            # Second Query
            FROM=make_query(
                # gaia_mags=True,
                user_cols=querystr['l1cols'], defaults='full',
                # Innermost Query
                FROM=make_query(
                    # gaia_mags=True,
                    defaults='full', inmostquery=True,
                    user_cols=querystr['l0cols'],
                    # random index for testing
                    random_index=_random_index,
    ))))

    if not as_qtable:
        df = Table(df)

    # -------------------------------
    # Saving

    if save_window:
        if isinstance(save_window, str):
            fpath = save_window
        else:
            fpath = 'output/window.fit'

        dump(df, fpath)

    if save_frame:
        if isinstance(save_frame, str):
            fpath = save_frame
        else:
            fpath = 'output/MasterFrame.pkl'

        dump([Frame, coordtransforms], fpath)


    # -------------------------------
    # Plotting
    plt.scatter(df['ra'], df['dec'], s=.1, c='b')
    plt.plot(orbit.icrs.ra, orbit.icrs.dec, c='k')
    plt.scatter(point.icrs.ra, point.icrs.dec, s=100, c='r')
    plt.set(title=f'Query around point {point.icrs.ra:.2},{point.icrs.dec:.2}',
            xlabel='\alpha [deg]', ylabel='\delta [deg]',
            savefig='plots/query_around_point{point.icrs.ra:.2}-{point.icrs.dec:.2}.png')

    # -------------------------------
    # Returning

    if return_transforms:
        return df, Frame, coordtransforms
    else:
        return df
# /def


##############################################################################
### query_along_orbit

def query_along_orbit(orbit, frame_point, frame_philims,
                      # points
                      points_inds=[], points_philims=None,
                      # options
                      frame_dt=1.9 * u.Myr, points_dt=None,
                      max_separation=.5 * u.kpc, return_transforms=True,
                      # for make_query
                      local=False,
                      _random_index=None,
                      #
                      distance=23.2 * u.kpc,
                      duststr='output/ps1dust_{}.dat',
                      # saving and logging
                      pprint=False,
                      save_window=False, save_frame=False,
                      logger=_LOGFILE, verbose=None):
    """

    Parameters
    ----------
    orbit : SkyCoord
        icrs orbit

    frame_point : SkyCoord

    points_philims : list
        list (as long points_inds) of philim lists (len 4)

    save_window : bool or str
        False : nothing
        True : save to 'output/window.fits'
        str : save at this location

    save_frame : bool or str
        save the main Frame and coordtransforms
        False : nothing
        True : save to 'output/window.fits'
        str : save at this location

    """
    # -------------------------------------------------------------------------
    # Preparing

    ## checks
    assert isinstance(save_window, (bool, str))
    assert isinstance(save_frame, (bool, str))

    ## points_philims
    if points_philims is None:
        points_philims = (frame_philims, ) * len(points_inds)
    elif len(points_philims) != len(points_inds):
        raise ValueError('points_philims do not match points_inds in length')

    ## dt
    if points_dt is None:
        points_dt = (frame_dt, ) * len(points_inds)
    elif len(points_dt) != len(points_inds):
        return ValueError('points_dt do not match points_inds in length')

    # -------------------------------------------------------------------------
    # preplotting

    plt.plot(orbit.icrs.ra, orbit.icrs.dec, c='k')
    plt.scatter(frame_point.icrs.ra, frame_point.icrs.dec, s=100, c='r')
    for i in points_inds:
        plt.scatter(orbit.icrs.ra[i], orbit.icrs.dec[i], s=50, c='k')
    plt.set(title='Query Along Orbit',
            xlabel='\alpha [deg]', ylabel='\delta [deg]',
            savefig='plots/query_along_point.png')

    # -------------------------------------------------------------------------
    # the master custom frame

    logfile.newsection(title=f'Starting Query:', div='=')

    df, MasterFrame, crd_xfm = query_at_point(
        orbit, frame_point, frame_philims, dt=frame_dt, is_nbody=False,
        max_separation=max_separation, return_transforms=True,
        as_qtable=False,
        local=local, _random_index=_random_index,
        pprint=pprint,
        save_window=False, save_frame=save_frame,
        logger=logfile, verbose=verbose
    )

    # -------------------------------------------------------------------------
    # Querying at points

    for i, (ind, pdt, philims) in enumerate(zip(points_inds, points_dt, points_philims)):

        point = orbit[ind]

        logfile.newsection(title=f'Query @ Point {i+1}:', div='=')

        pointdf = query_at_point(
            orbit, point, philims, dt=pdt, is_nbody=False,
            max_separation=max_separation, return_transforms=False,
            as_qtable=False,
            local=local, _random_index=_random_index,
            pprint=pprint,
            logger=logfile, verbose=verbose)

        # xmatch the catalogs to identify duplicates
        _, idx2, _ = xmatch(
            df, pointdf,
            colRA1='ra', colDec1='dec', epoch1=2015.5,
            colRA2='ra', colDec2='dec', epoch2=2015.5)

        # finding distinct values in pointdf
        inds = np.arange(len(pointdf))         # full index array for comparison
        nm = list(set(inds).difference(idx2))  # not matched values

        # merging tables
        df = vstack([df, pointdf[nm]])

    # -------------------------------------------------------------------------
    # modifying table

    register_stream_frame(crd_xfm.R_icrs_to_cust)

    # replace the custom coordinates by values from MasterFrame
    sc = SkyCoord(ra=df['ra'], dec=df['dec'],
                  pm_ra_cosdec=df['pmra'], pm_dec=df['pmdec'],
                  frame='icrs', representation_type='spherical'
                  ).transform_to(MasterFrame)
    sc.phi1.wrap_angle = '180d'

    df['phi1'] = sc.phi1.to(deg)
    df['phi2'] = sc.phi2.to(deg)
    df['pmphi1'] = sc.pm_phi1_cosphi2
    df['pmphi2'] = sc.pm_phi2

    df = QTable(df)

    df['prlx'][df['prlx'] < 0] = np.NaN

    table_utils.add_color_col(df, 'g', 'r')
    table_utils.add_color_col(df, 'g', 'i')
    table_utils.add_abs_pm_col(df, 'pmra', 'pmdec')

    logfile.write("adding color columns", "did g-r", "did g-i", "did |pm|",
                  sep='\n\t', end='\n')

    # -------------------------------------------------------------------------
    # Dust

    # distance = pal5['sc'][0].distance.to_value(u.kpc)

    ps1dust = load_dust_gri(duststr, df=df, distance=distance)

    ## Making color corrections
    logfile.write("Dextincting",
                   "g -> g dx", "r -> r dx", "i -> i dx",
                   "g-r -> g-r dx", "g-i -> g-i dx", "r-i -> r-i dx",
                   sep='\n\t', end='\n')

    df['g dx'] = df['g'] - ps1dust['g']
    df['r dx'] = df['r'] - ps1dust['r']
    df['i dx'] = df['i'] - ps1dust['i']

    df['g-r dx'] = df['g dx'] - df['r dx']
    df['g-i dx'] = df['g dx'] - df['i dx']
    df['r-i dx'] = df['r dx'] - df['i dx']

    for c in ('g', 'r', 'i', 'g-r', 'g-i', 'r-i'):
        setattr(df[c + ' dx'].info, 'description', 'de-extincted w/ mwdust')

    ## Adding in Conversion to MegaCam
    df['g MC'] = MCG1_PS1.G_MP9401(df, g='g dx', r='r dx', gmi='g-r dx')
    df['r MC'] = MCG1_PS1.R_MP9601(df, g='g dx', r='r dx', gmi='g-r dx')
    df['g-r MC'] = df['g MC'] - df['r MC']

    for c in ('g MC', 'r MC', 'g-r MC'):
        setattr(df[c].info, 'description', 'de-extincted & converted to MegaCam')

    # -------------------------------------------------------------------------
    # Saving

    df = QTable(df)

    if save_window:  # True or str
        if isinstance(save_window, str):
            fpath = save_window
        else:
            fpath = 'output/window.fits'
        logfile.write(f'saving to {fpath}')
        df.write(fpath, format='fits', overwrite=True)

    # -------------------------------------------------------------------------
    # postplotting

    plt.scatter(df['ra'].to_value(u.deg), df['dec'].to_value(u.deg), s=.1, c='k')
    plt.plot(orbit.icrs.ra, orbit.icrs.dec, c='k')
    plt.scatter(frame_point.icrs.ra, frame_point.icrs.dec, s=100, c='r')
    for i in points_inds:
        plt.scatter(orbit.icrs.ra[i], orbit.icrs.dec[i], s=50, c='k')
    plt.set(title='Query Along Orbit',
            xlabel='\alpha [deg]', ylabel='\delta [deg]',
            savefig='plots/query_along_point.png')


    # plt.scatter(df['phi1'], df['phi2'], s=.1, c='k', fig='new')
    # plt.plot(orbit.icrs.ra, orbit.icrs.dec, c='k')
    # plt.scatter(frame_point.icrs.ra, frame_point.icrs.dec, s=100, c='r')
    # for i in points_inds:
    #     plt.scatter(orbit.icrs.ra[i], orbit.icrs.dec[i], s=50, c='k')
    # plt.set(title='Query Along Orbit',
    #         xlabel='\alpha [deg]', ylabel='\delta [deg]',
    #         savefig='plots/query_along_point.png')

    # -------------------------------------------------------------------------
    # returning

    return df, MasterFrame, crd_xfm
# /def


# ##############################################################################
# ### DONE
