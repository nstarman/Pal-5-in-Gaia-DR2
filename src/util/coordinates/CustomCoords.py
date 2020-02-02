# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : Custom Coordinates
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Custom Coordinates.

TODO
----
replace custom_to_radec with version from galpy.util.bovy_coords
    and delete custom_to_radec, from mybovy_coords

"""

__author__ = "Nathaniel Starkman"

#############################################################################
# IMPORTS

# GENERAL
from numpy.linalg import norm

# astropy
from astropy import units as u, constants as consts
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, frame_transform_graph
from astropy.coordinates.representation import CartesianDifferential
from astropy.coordinates.matrix_utilities import (
    rotation_matrix,
    matrix_product,
    matrix_transpose,
)

# galpy
from galpy.util.bovy_coords import radec_to_custom, custom_to_radec

# CUSTOM
from astroPHD import ObjDict, LogFile

# PROJECT-SPECIFIC
from ..stats import mean as average


#############################################################################
# PARAMETERS

_LOGFILE = LogFile(header=False)  # LogPrint, which is compatible with LogFile


#############################################################################
# CODE
#############################################################################

#############################################################################
# Rotated Coordinate Frame Class


class StreamFrame(coord.BaseCoordinateFrame):
    """StreamFrame.

    A Heliocentric spherical coordinate system defined to linearize
    a stream about a point, using the angular momentum at that point.

    http://docs.astropy.org/en/stable/generated/examples/coordinates/
        plot_sgr-coordinate-frame.html
        #sphx-glr-generated-examples-coordinates-plot-sgr-coordinate-frame-py

    """

    default_representation = coord.SphericalRepresentation
    default_differential = coord.SphericalCosLatDifferential

    frame_specific_representation_info = {
        coord.SphericalRepresentation: [
            coord.RepresentationMapping("lon", "phi1"),
            coord.RepresentationMapping("lat", "phi2"),
            coord.RepresentationMapping("distance", "distance"),
        ],
        coord.SphericalCosLatDifferential: [
            coord.RepresentationMapping("d_lon_coslat", "pm_phi1_cosphi2"),
            coord.RepresentationMapping("d_lat", "pm_phi2"),
            coord.RepresentationMapping("d_distance", "radial_velocity"),
        ],
        coord.SphericalDifferential: [
            coord.RepresentationMapping("d_lon", "pm_phi1"),
            coord.RepresentationMapping("d_lat", "pm_phi2"),
            coord.RepresentationMapping("d_distance", "radial_velocity"),
        ],
    }

    frame_specific_representation_info[
        coord.UnitSphericalRepresentation
    ] = frame_specific_representation_info[coord.SphericalRepresentation]
    frame_specific_representation_info[
        coord.UnitSphericalCosLatDifferential
    ] = frame_specific_representation_info[coord.SphericalCosLatDifferential]
    frame_specific_representation_info[
        coord.UnitSphericalDifferential
    ] = frame_specific_representation_info[coord.SphericalDifferential]


# /class


##############################################################################


def register_stream_frame(R_icrs_to_cust, logger=_LOGFILE, verbose=None):
    """Register StreamCoords into the Astropy frame_transform_graph.

    Parameters
    ----------
    R_icrs_to_cust: rotation matrix
        rotates from Cartesian ICRS coords to Cartesian Custom Coords

    """
    logger.report("registered StreamCoordFrame", verbose=verbose)

    @frame_transform_graph.transform(
        coord.StaticMatrixTransform, coord.ICRS, StreamFrame
    )
    def icrs_to_cust():
        """Register matrix transformation.

        Does ICRS spherical to ICRS-rotated custom coordinates

        """
        return R_icrs_to_cust

    # /def

    @frame_transform_graph.transform(
        coord.StaticMatrixTransform, StreamFrame, coord.ICRS
    )
    def cust_to_icrs():
        """Register matrix transformation.

        Does ICRS-rotated custom to ICRS spherical coordinates

        """
        return matrix_transpose(R_icrs_to_cust)

    # /def

    return

# /def


# -----------------------------------------------------------------------------


def reference_to_skyoffset_matrix(skyoffset_frame):
    """Make rotation matrix from reference frame and a skyoffset frame.

    this function is useful after using astropy to make a rotation matrix
    *modified from reference_to_skyoffset function in astropy

    Parameters
    ----------
    skyoffset_frame : SkyCoord frame
        the rotated frame
        note: SkyCoords are compatible

    Returns
    -------
    R : numpy array
        a 3d rotation matrix

    """
    origin = skyoffset_frame.origin.spherical

    mat1 = rotation_matrix(-skyoffset_frame.rotation, "x")
    mat2 = rotation_matrix(-origin.lat, "y")
    mat3 = rotation_matrix(origin.lon, "z")

    R = matrix_product(mat1, mat2, mat3)
    return R


# /def


# -----------------------------------------------------------------------------


def make_stream_frame_transform(
    orbit,
    point,
    method="path",
    dt=2 * u.Myr,
    is_nbody=False,
    max_separation=0.5 * u.kpc,
    logger=_LOGFILE,
    verbose=None,
):
    """Make_stream_frame_transform.

    Parameters
    ----------
    orbit: the orbits
        columns are [x, y, z]
    point: SkyCoord
        starting point
        ex: orbit.icrs.[0]
    method : str
        options:
            'path' : using the path of the orbit  (TODO average path of stream)

    dt : time
        only used if *method* = 'path'

    """
    logger.report("making make_stream_frame_transform", verbose=verbose)

    if method == "path":

        if is_nbody:  # need to redifine point as average
            idx = orbit.separation_3d(point) < max_separation
            point = average(orbit[idx])

        # offset from point
        motion_coords = point.apply_space_motion(dt=dt)

        # snap to actual orbit
        if is_nbody:  # an average of nearby things
            idx = orbit.separation_3d(motion_coords) < max_separation
            offset_coords = average(orbit[idx])
        else:
            idx, sep2d, dist3d = motion_coords.match_to_catalog_sky(orbit)
            offset_coords = orbit[idx]

        # get angle of rotation of rotated frame
        rotation = point.position_angle(offset_coords) - 90 * u.deg

        # the custom frame
        # made by astropy and then used to extract the correct rotation matrix
        custframe = point.skyoffset_frame(rotation=rotation)

        # icrs to custom rotation matrix
        R_icrs_to_cust = reference_to_skyoffset_matrix(custframe)

    else:
        raise ValueError()

    logger.report(R_icrs_to_cust)

    # NGP angles
    theta_ngp, _ = (
        radec_to_custom(180.0, 90.0, T=R_icrs_to_cust, degree=True) * u.deg
    )
    ra_ngp, dec_ngp = (
        custom_to_radec(180.0, 90.0, T=R_icrs_to_cust, degree=True) * u.deg
    )

    coordtransforms = ObjDict(
        "coordinate transforms",
        # Rcust_gal=Rcust_gal,
        R_cust_eq=R_icrs_to_cust,
        R_icrs_to_cust=R_icrs_to_cust,
        # Lvec=Lvec,
        # nrpnt=nrpnt,
        ra_ngp=ra_ngp,
        dec_ngp=dec_ngp,
        theta_ngp=theta_ngp,
        point=point,
    )

    return coordtransforms


# /def


# -----------------------------------------------------------------------------


def make_stream_frame(
    orbit,
    point,
    method="path",
    dt=2 * u.Myr,
    is_nbody=False,
    max_separation=0.5 * u.kpc,
    return_transforms=False,
    logger=_LOGFILE,
    verbose=None,
):
    """Makes and registers StreamCoord reference frame.

    Parameters
    ----------
    orbit :
        the orbit to use to rotate
    point :
        point in cartesian
        starting point = orbit['gal'].cartesian.xyz.T[0]
    method : str
        options:
            'angular momentum' : the angular momentum  (TODO re-add support)
            'path' : using the path of the orbit
    is_nbody : bool, optional
        default: False
        the point is the sp
        snap to nearest point and use averaging
        if it's not an orbit,

    Returns
    -------
    StreamFrame : astropy Coordinate Frame

    (if return_transforms is True)
    coordtransforms : ObjDict

    Exceptions
    ----------
    ValueError : if method is not in allowed options

    """
    if method not in ("path", "angular momentum"):
        raise ValueError("method not supported. See documentation.")

    coordtransforms = make_stream_frame_transform(
        orbit,
        point,
        method=method,
        dt=dt,
        is_nbody=is_nbody,
        max_separation=max_separation,
        logger=logger,
        verbose=verbose,
    )

    # register frame into astropy's frame_transform_graph
    register_stream_frame(
        coordtransforms.R_icrs_to_cust, logger=logger, verbose=verbose
    )

    if return_transforms:
        return StreamFrame, coordtransforms
    return StreamFrame


# /def


##############################################################################
# END
