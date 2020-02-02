#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#
# TITLE   : frameattrs
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# -----------------------------------------------------------------------------

# Docstring and Metadata
"""Frame Attributes
"""

__author__ = "Nathaniel Starkman"

#############################################################################
# IMPORTS

# General
import warnings
import numpy as np

# astropy
from astropy import units as u
from astropy import coordinates as coord
from astropy.coordinates.baseframe import FrameMeta, _get_repr_cls
from astropy.coordinates.representation import (
    MetaBaseRepresentation,
    MetaBaseDifferential,
)


#############################################################################
# CODE
##############################################################################


def frameattrs(
    frame, representation=None, differential=None, notrecog2custom=False
):
    """frameattrs.

    frame: str or SkyCoord
    representation: str or astropy representation class
    differential: str or astropy differential representation class

    TODO use sc.get_representation_component_names() instead of a dictionary
    TODO not always kpc & deg!!!

    """
    reps = {
        coord.CartesianRepresentation: "cartesian",
        coord.SphericalRepresentation: "spherical",
        coord.CylindricalRepresentation: "cylindrical",
        coord.RadialRepresentation: "radial",
    }

    diffs = {
        coord.SphericalCosLatDifferential: "sphericalcoslat",
        coord.SphericalDifferential: "spherical",
        coord.CartesianDifferential: "cartesian",
        coord.CylindricalDifferential: "cylindrical",
        coord.RadialDifferential: "radial",
    }
    recogframe = ("icrs", "galactic", "galactocentric")  # known frame

    # +--------------- Frame ---------------+
    if isinstance(frame, str):
        framename = frame
        if framename not in (*recogframe, "custom"):
            raise ValueError(f"{framename} is not a recognized frame")

    # elif isinstance(frame, FrameMeta):  # TODO how detect custom frame?
    else:
        framename = frame.name

        if framename not in recogframe:
            if not notrecog2custom:
                raise ValueError(f"{framename} is not a recognized frame")
            framename = "custom"
    # else:  # TODO how to detect custom frames?
    #     raise ValueError(f'frame is not <str> or Frame')

    # +--------------- Representation ---------------+
    if representation is None:
        # if not isinstance(frame, FrameMeta):  # TODO how detect custom frame?
        #     raise ValueError('representation cannot be None if frame is str')

        if frame.representation_type not in reps.keys():
            raise ValueError(
                (
                    f"{frame.representation_type} not representation of"
                    f"{reps.values()}"
                )
            )

        representation = reps[frame.representation_type]

    elif isinstance(representation, str):
        if representation not in reps.values():
            raise ValueError(f"representation not in {reps.values()}")
        pass  # already a str

    elif isinstance(representation, MetaBaseRepresentation):
        if representation not in reps.keys():
            raise ValueError(
                f"{representation} not representation of {reps.values()}"
            )

        representation = reps[representation]

    # +--------------- Differential ---------------+
    if differential is None:
        # if not isinstance(frame, FrameMeta):  # TODO how detect custom frame?
        #     raise ValueError('differential cannot be None if frame is str')

        if frame.differential_type not in diffs.keys():
            raise ValueError(
                (
                    f"{frame.differential_type} not representation of"
                    f"{diffs.values()}"
                )
            )

        differential = diffs[frame.differential_type]

    if isinstance(differential, str):
        if differential not in diffs.values():
            raise ValueError(f"differential not in {diffs.values()}")
        pass

    elif isinstance(differential, MetaBaseDifferential):
        if differential not in diffs.keys():
            raise ValueError(
                f"{differential} not representation of {diffs.values()}"
            )

        differential = diffs[differential]

    # +--------------- Attributes & Units ---------------+
    attrs = []
    units = {}

    # ICRS
    if framename.lower() == "icrs":
        # Representation
        if representation == "cartesian":
            attrs.extend(("x", "y", "z"))
            units.update({"x": u.kpc, "y": u.kpc, "z": u.kpc})
        elif representation == "spherical":
            attrs.extend(("ra", "dec", "distance"))
            units.update({"ra": u.deg, "dec": u.deg, "distance": u.kpc})
        elif representation == "cylindrical":
            attrs.extend(("rho", "phi", "z"))
            units.update({"rho": u.kpc, "phi": u.deg, "z": u.kpc})
        elif representation == "radial":
            raise Exception

        # Differential
        if differential == "cartesian":
            attrs.extend(("v_x", "v_y", "v_z"))
            units.update(
                {"v_x": u.km / u.s, "v_y": u.km / u.s, "v_z": u.km / u.s}
            )
        elif differential == "spherical":
            attrs.extend(("pm_ra", "pm_dec", "radial_velocity"))
            units.update(
                {
                    "pm_ra": u.mas / u.yr,
                    "pm_dec": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "sphericalcoslat":
            attrs.extend(("pm_ra_cosdec", "pm_dec", "radial_velocity"))
            units.update(
                {
                    "pm_ra_cosdec": u.mas / u.yr,
                    "pm_dec": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "cylindrical":
            attrs.extend(("d_rho", "d_phi", "d_z"))
            units.update(
                {"d_rho": u.km / u.s, "d_phi": u.mas / u.yr, "d_z": u.km / u.s}
            )
        elif differential == "radial":
            raise Exception

    # Galactic
    elif framename.lower() == "galactic":
        # Representation
        if representation == "cartesian":
            attrs.extend(("u", "v", "w"))
            units.update({"u": u.kpc, "v": u.kpc, "w": u.kpc})
        elif representation == "spherical":
            attrs.extend(("l", "b", "distance"))
            units.update({"l": u.deg, "b": u.deg, "distance": u.kpc})
        elif representation == "cylindrical":
            attrs.extend(("rho", "phi", "z"))
            units.update({"rho": u.kpc, "phi": u.deg, "z": u.kpc})
        elif representation == "radial":
            raise Exception

        # Differential
        if differential == "cartesian":
            attrs.extend(("U", "V", "W"))
            units.update({"U": u.km / u.s, "V": u.km / u.s, "W": u.km / u.s})
        elif differential == "spherical":
            attrs.extend(("pm_l", "pm_b", "radial_velocity"))
            units.update(
                {
                    "pm_l": u.mas / u.yr,
                    "pm_b": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "sphericalcoslat":
            attrs.extend(("pm_l_cosb", "pm_b", "radial_velocity"))
            units.update(
                {
                    "pm_l_cosb": u.mas / u.yr,
                    "pm_b": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "cylindrical":
            attrs.extend(("d_rho", "d_phi", "d_z"))
            units.update(
                {"d_rho": u.km / u.s, "d_phi": u.mas / u.yr, "d_z": u.km / u.s}
            )
        elif differential == "radial":
            raise Exception

    # Galactocentric
    elif framename.lower() == "galactocentric":
        # Representation
        if representation == "cartesian":
            attrs.extend(("x", "y", "z"))
            units.update({"x": u.kpc, "y": u.kpc, "z": u.kpc})
        elif representation == "spherical":
            attrs.extend(("lon", "lat", "distance"))
            units.update({"lon": u.deg, "lat": u.deg, "distance": u.kpc})
        elif representation == "cylindrical":
            attrs.extend(("rho", "phi", "z"))
            units.update({"rho": u.kpc, "phi": u.deg, "z": u.kpc})
        elif representation == "radial":
            raise Exception

        # Differential
        if differential == "cartesian":
            attrs.extend(("v_x", "v_y", "v_z"))
            units.update(
                {"v_x": u.km / u.s, "v_y": u.km / u.s, "v_z": u.km / u.s}
            )
        elif differential == "spherical":
            attrs.extend(("pm_lon", "pm_lat", "radial_velocity"))
            units.update(
                {
                    "pm_lon": u.mas / u.yr,
                    "pm_lat": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "sphericalcoslat":
            attrs.extend(("pm_lon_coslat", "pm_lat", "radial_velocity"))
            units.update(
                {
                    "pm_lon_coslat": u.mas / u.yr,
                    "pm_lat": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "cylindrical":
            attrs.extend(("d_rho", "d_phi", "d_z"))
            units.update(
                {"d_rho": u.km / u.s, "d_phi": u.mas / u.yr, "d_z": u.km / u.s}
            )
        elif differential == "radial":
            raise Exception

    # Custom
    elif framename.lower() == "custom":
        # Representation
        if representation == "cartesian":
            raise Exception
        elif representation == "spherical":
            attrs.extend(("phi1", "phi2", "distance"))
            units.update({"phi1": u.deg, "phi2": u.deg, "distance": u.kpc})
        elif representation == "cylindrical":
            raise Exception
        elif representation == "radial":
            raise Exception

        # Differential
        if differential == "cartesian":
            raise Exception
        elif differential == "spherical":
            attrs.extend(("pm_phi1", "pm_phi2", "radial_velocity"))
            units.update(
                {
                    "pm_phi1": u.mas / u.yr,
                    "pm_phi2": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "sphericalcoslat":
            attrs.extend(("pm_phi1_cosphi2", "pm_phi2", "radial_velocity"))
            units.update(
                {
                    "pm_phi1_cosphi2": u.mas / u.yr,
                    "pm_phi2": u.mas / u.yr,
                    "radial_velocity": u.km / u.s,
                }
            )
        elif differential == "cylindrical":
            raise Exception
        elif differential == "radial":
            raise Exception

    else:
        raise ValueError(
            f"{frame} not in (ICRS, Galactic, Galactocentric, Custom)"
        )

    return attrs, units


# /def

# ----------------------------------------------------------------------------


def transform_sc(
    sc,
    frame=None,
    representation=None,
    representation_type=None,
    differential=None,
    differential_type=None,
    wrap_angle=None,
    **wrap_angles,
):
    """Transform a SkyCoord to the given frame, representation/differential_type.

    Arguments
    ---------
    sc: SkyCoord
    frame: str, None, astropy frame   (default: None)
        the desired frame.
        Can accept any astropy-compatible frame on the transform graph
        None: do not transform
        str: must be registered on the transform graph
            standard:
                icrs, galactic, galactocentric
    representation: None, str or astropy representation class
        None: get from
        available str:
            spherical, cartesian, cylindrical, radial
    differential: str or astropy differential representation class
        available str:
            sphericalcoslat, spherical, cartesian, cylindrical, radial
    wrap_angle: list
        (name, wrap_angle) or [(name, wrap_angle), (name, wrap_angle), ...]
        superseded by kwargs
        __deprecated__

    Other Parameters
    ----------------
    angle and wrapping

    """
    recogframe = ("icrs", "galactic", "galactocentric")  # known frame

    # +--------------- Frame ---------------+
    _transform = True  # whether to transform

    # default frame
    if frame is None:
        frame = sc.frame
        _transform = False  # already right frame

    # now checking valid frame
    if isinstance(frame, str):
        if frame == sc.frame.name:
            _transform = False
            pass
        elif frame not in recogframe:
            raise ValueError(f"frame is not recognized")
    elif isinstance(frame, FrameMeta):
        pass
    else:  # TODO how to detect custom frames?
        # raise TypeError(f'frame is not <str> or Frame')
        pass

    # +--------------- Representation ---------------+
    if representation_type is not None:
        representation = representation_type

    if representation is None:
        if isinstance(frame, str):
            warnings.warn("no representation_type specified, keeping current")
            representation_type = None
        else:
            try:
                _get_repr_cls(frame.representation_type)
            except ValueError:  # frame is a baseframe, keep original
                representation_type = sc.representation_type
            else:
                representation_type = frame.representation_type

    elif isinstance(representation, str):
        if representation == "spherical":
            representation_type = coord.SphericalRepresentation
        elif representation == "cartesian":
            representation_type = coord.CartesianRepresentation
        elif representation == "cylindrical":
            representation_type = coord.CylindricalRepresentation
        elif representation == "radial":
            representation_type = coord.RadialRepresentation
        else:
            raise ValueError

    elif isinstance(representation, MetaBaseRepresentation):
        representation_type = representation

    else:
        raise TypeError

    # +--------------- Differential ---------------+
    if differential_type is not None:
        differential = differential_type

    if differential is None:  # no differential
        if isinstance(frame, str):  # keep original differential
            warnings.warn("no differential_type specified, keeping current")
            differential_type = None
        else:  # transform to frame differential
            try:
                _get_repr_cls(frame.differential_type)
            except ValueError:  # frame is a baseframe, keep original
                differential_type = sc.differential_type
            else:
                differential_type = frame.differential_type

    elif isinstance(differential, str):  # override frame differential
        if differential == "sphericalcoslat":
            differential_type = coord.SphericalCosLatDifferential
        elif differential == "spherical":
            differential_type = coord.SphericalDifferential
        elif differential == "cartesian":
            differential_type = coord.CartesianDifferential
        elif differential == "cylindrical":
            differential_type = coord.CylindricalDifferential
        elif differential == "radial":
            differential_type = coord.RadialDifferential
        else:
            raise ValueError

    elif isinstance(differential, MetaBaseDifferential):
        differential_type = differential

    else:
        raise TypeError

    # +--------------- Converting ---------------+
    if _transform:
        sc = sc.transform_to(frame)
    if representation_type is not None:
        sc.representation_type = representation_type
    if differential_type is not None:
        sc.differential_type = differential_type

    # TODO deprecate
    if wrap_angle is not None:
        if np.isscalar(wrap_angle[0]):  # just (name, wrapangle)
            setattr(getattr(sc, wrap_angle[0]), "wrap_angle", wrap_angle[1])
        else:  # list [(name, wrapangle), (name, wrapangle), ...]
            for name, angle in wrap_angle:
                setattr(getattr(sc, name), "wrap_angle", angle)

    for name, angle in wrap_angles.items():
        setattr(getattr(sc, name), "wrap_angle", angle)

    return sc


# /def

# ----------------------------------------------------------------------------


def convert_frame_and_repr(
    sc, frame=None, representation=None, differential=None
):
    """convert_frame_and_repr."""
    return transform_sc(
        sc,
        frame=frame,
        representation=representation,
        differential=differential,
    )


convert_frame_and_repr.__doc__ = transform_sc.__doc__
# /def

# ----------------------------------------------------------------------------


def convert_repr(sc, representation=None, differential=None):
    """convert_repr."""
    return transform_sc(
        sc,
        frame=None,
        representation=representation,
        differential=differential,
    )


# /def


##############################################################################
# END
