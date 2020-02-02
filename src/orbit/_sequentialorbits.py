# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : Sequential Orbits
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""SequentialOrbits."""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# GENERAL
import numpy as np
from typing import Optional, Union

# Galpy
from galpy.potential import MWPotential2014, Potential
from galpy.orbit import Orbit

# Astropy
from astropy.coordinates import concatenate, SkyCoord
from astropy import units as u


# CUSTOM
from astroPHD.community import starkplot as plt
from astroPHD.community.starkplot import mpl_decorator


# PROJECT-SPECIFIC
from .OrbitKeeper import OrbitKeeper

from .exceptions import warn, integrationWarning
from ..util import astrarray
from ..util.pickle import dump as _dump, load as _load


###############################################################################
# Code
###############################################################################


class SequentialOrbits(object):
    """SequentialOrbit.

    A class which holds many segments of one orbit.

    """

    def __new__(cls, *args, **kw):
        """Make new SequentialOrbits.

        used for copying SequentialOrbits

        """
        # assume first argument is a real argument
        if len(args) > 0:  # vxvv first argument
            vxvv = args[0]
        else:  # vxvv kwarg
            vxvv = kw.get("vxvv", None)

        # # for making a blank instance  # TODO obviate
        # if vxvv is None:  # needed for pickling, but calls __new__ twice
        #     return super().__new__(cls)  # empty new instance

        # for copying
        if isinstance(vxvv, SequentialOrbits):  # copies by __new__() again :(
            instance = cls.__new__(cls)
            instance.__init__(
                vxvv=vxvv._data,  # using internal OrbitKeeper
                pot=vxvv._pot,
                ro=vxvv._ro,
                vo=vxvv._vo,
                zo=vxvv._zo,
                solarmotion=vxvv._solarmotion,
            )
            return instance

        else:
            self = super().__new__(cls)  # make blank instnce to start __init__
        return self

    # /def

    def __init__(
        self,
        vxvv=None,
        lb: bool = False,
        radec: bool = False,
        ro: Optional[u.Quantity] = None,
        vo: Optional[u.Quantity] = None,
        zo: Optional[u.Quantity] = None,
        solarmotion: Union[str, list] = "stream",
        tstart: Optional[u.Quantity] = None,
        pot: Optional[Potential] = None,
    ):
        """Instantiate SequentialOrbits.

        Parameters
        ----------
        vxvv: Orbits or any Orbits-compatible input
            the initial conditions
        lb:
        radec:
        ro: Quantity  (default = 8 kpc)
            distance from vantage point to GC (kpc)
        vo: Quantity  (default = 220 km/s)
            circular velocity at ro (km/s)
        zo: Quantity  (default = 25 pc)
            offset toward the NGP of the Sun wrt the plane (kpc)
        solarmotion:  str or array
            str: 'stream', 'local', 'hogg' or 'dehnen', or 'schoenrich',
            array: value in [-U,V,W]. can be Quantity
        tstart:  Quantity  (default = 0 gyr)
            start time of first orbit segment
        pot:  Galpy potential
            potential

        TODO
        ------
        allow instantiation from IndexDict

        """
        # +---------- Galactic Properties ----------+
        self._ro = (
            ro if ro is not None else 8.0 * u.kpc
        )  # TODO set by default?
        self._vo = vo if vo is not None else 220.0 * u.km / u.s
        self._zo = zo if zo is not None else 25 * u.pc
        # solarmotion
        if isinstance(solarmotion, str):  # TODO support all galpy defaults
            if solarmotion.lower() == "stream":
                self._solarmotion = [-11.1, 24.0, 7.25] * u.km / u.s
            elif solarmotion.lower() == "local":
                self._solarmotion = [-11.1, 12.24, 7.25] * u.km / u.s
        else:
            self._solarmotion = solarmotion

        # +---------- Orbit Properties ----------+
        # starting time
        t0: u.Quantity = 0 * u.Gyr if tstart is None else tstart

        # making orbits
        if isinstance(vxvv, OrbitKeeper):  # internal use only
            self._data = vxvv
        # elif isinstance(vxvv, IndexDict):  # TODO
        #     self._data = OrbitKeeper(None, vxvv)
        #     self._data._ind = vxvv._ind
        #     self._data.direction = vxvv.direction
        #     self._data.t0 = vxvv[vxvv.minkey].time()[0]  # overriding t0
        #     self._data.tref = vxvv[vxvv.maxkey].time()[0]  # overriding t0
        #     # self._data._times =
        #     # self._data._bounds =
        elif isinstance(vxvv, Orbit):
            self._data = OrbitKeeper(vxvv, t0)
        else:  # accepts: lists, SkyCoord
            orbit = Orbit(
                vxvv=vxvv,
                ro=self._ro,
                vo=self._vo,
                zo=self._zo,
                solarmotion=self._solarmotion,
                radec=radec,
                uvw=False,
                lb=lb,
            )
            self._data = OrbitKeeper(orbit, t0)

        # potential
        # defaults to MWPotential2014
        self._pot = pot if pot is not None else MWPotential2014

        return

    # /def

    # +---------- time ----------+
    @property
    def times(self):
        """Return all times."""
        return self._data._times

    @property
    def time(self):
        """Return current time."""
        return self._data.time

    @time.setter
    def time(self, value):
        """Return current time."""
        self._data.time = value

    @property
    def t(self):
        """__deprecated__ current time."""
        return self._data.time

    def set_time(self, time):
        """Set current time array."""
        self.time = time

    # +---------- bounds ----------+
    # bounds are kept mostly hidden
    @property
    def _allbounds(self):
        """All bounds."""
        return self._data.allbounds

    @property
    def _bounds(self):
        """Return current bounds."""
        return self._data.bounds

    @_bounds.setter
    def _bounds(self, value):
        """Return current bounds."""
        self._data.bounds = value

    def _set_bounds(self, bounds):
        """Set current bounds."""
        self._data.bounds = bounds

    # +---- Orbits  ----+
    @property
    def _orbits(self):
        """All orbits."""
        return self._data._orbits

    @property
    def orbit(self):
        """Return current orbit."""
        return self._data.orbit

    @property
    def o(self):
        """Return current orbit."""
        return self._data.orbit

    @property
    def _vxvv(self):
        """Return current orbit initial phase parameters."""
        return self._data.orbit.vxvv

    @property
    def _init_vxvv(self):
        """Return initial phase parameters of 1st parameters."""
        return self._data._orbits[0].vxvv

    @property
    def segmentlist(self):
        """Orbit segments list."""
        return self._data.keyslist()

    @property
    def numsegs(self):
        """Return number of orbit segments."""
        # can't use self._data._ind b/c deleted orbits
        return len(self._data.keyslist())

    @property
    def numorbs(self):
        """Return number of orbits (per segment)."""
        return len(self._data.orbit)

    # +---------- potential ----------+
    @property
    def potential(self):
        """Potential."""
        return self._pot

    @potential.setter
    def potential(self, value):
        self._pot = value

    ##############################

    def __getitem__(self, ind):
        """__getitem__.

        slices operate on self.segmentlist!
        TODO support reindex option so have 0, ...

        """
        def _helper(ind0):
            if ind0 in (Ellipsis, None):
                return self._data
            elif isinstance(ind0, slice):
                print("slice", self.segmentlist[ind0])
                return self._data[self.segmentlist[ind0]]
            else:  # many indices return SequentialOrbits
                return self._data[ind0]

        # +----- single dimensional -----+
        if (
            np.isscalar(ind)
            or isinstance(ind, slice)
            or ind in (Ellipsis, None)
        ):
            if np.isscalar(ind):  # single index return Orbits
                return self._orbits[ind]  # list b/c IndexDict __getitem__
            else:
                vxvv = _helper(ind)
                orb = SequentialOrbits(
                    vxvv=vxvv,
                    ro=self._ro,
                    vo=self._vo,
                    zo=self._zo,
                    solarmotion=self._solarmotion,
                    pot=self.potential,
                )
                return orb

        # +----- multi-dimensional -----+
        elif len(ind) <= 3:
            # +----- single 0th index -----+
            if np.isscalar(ind[0]):  # single index return Orbits
                if len(ind) == 2:  # get Orbits subset
                    return self._orbits[ind]  # (tuple)
                elif len(ind) == 3:  # get Orbits subset at specific times
                    return self._orbits[ind[:2]].SkyCoord(ind[2])

            # +----- slice/array 0th index -----+
            else:
                vxvv = _helper(ind[0])

                if len(ind) == 2:  # TODO SequentialOrbits
                    vxvv = vxvv[:, ind[1]]  # (tuple)
                elif len(ind) == 3:
                    return vxvv[:, ind[1]].SkyCoord(ind[2])
                orb = SequentialOrbits(
                    vxvv=vxvv,
                    ro=self._ro,
                    vo=self._vo,
                    zo=self._zo,
                    solarmotion=self._solarmotion,
                    pot=self.potential,
                )
                return orb
        else:
            raise ValueError()

    # /def

    def __getattr__(self, name: str):
        """__getattr__.

        get or evaluate an attribute for this CompositeOrbit
        called if not already defined in Class

        Parameters
        ----------
        name: str
            name of the attribute

        Returns
        -------
        attr
            if the attribute is callable, a function to evaluate the attribute
            for each Orbit; otherwise a list of attributes

        """
        # Catch all plotting functions
        if "plot" in name:

            @mpl_decorator(fig=None)
            def _plot(*args, **kw):
                if len(args) == 0:
                    if "d1" not in kw:
                        kw["d1"] = "t"
                elif isinstance(args[0], str):
                    kw["d1"] = args[0]
                    args = args[1:]

                ind = self._data.keyslist()[0]

                # Hack to get the plot properties correct
                # since overplot doesn't make them
                fig, ax = plt.gcf(), plt.gca()
                plt.figure()
                getattr(self._orbits[ind], name)(*args, **kw)
                xl, yl = plt.gca().get_xlabel(), plt.gca().get_ylabel()
                plt.close()
                plt.figure(fig.number)
                plt.sca(ax)

                # actually plotting
                [
                    getattr(self._orbits[i], name)(*args, overplot=True, **kw)
                    for i in self._data.keys()
                ]

                # axisLabels(ax=None, x=xl, y=yl, units=False, **kw)
                ax.set_xlabel(xl)
                ax.set_ylabel(yl)

                return None

            # Assign documentation
            _plot.__doc__ = self._orbits[0][0].__getattribute__(name).__doc__
            return _plot

        # Do rest of functions

        # TODO don't make a new orbit every time
        attribute = getattr(Orbit(), name)

        if callable(attribute):

            def _func(t=None, *args, **kwargs):
                res = [
                    getattr(self._orbits[i], name)(t=time, *args, **kwargs)
                    for i, time in self.iterator(t)
                ]
                if len(res) == 1:
                    return res[0]
                else:
                    try:
                        return astrarray(res)
                    except AttributeError:
                        return res

            _func.__doc__ = attribute.__doc__
            return _func
        else:
            return [getattr(orbit, name) for orbit in self._orbits]

    # /def

    # +------------------------------+
    # Equality:
    def __eq__(self, other):
        """Test for equality.

        in: ro, vo, zo, solarmotion, potential, and orbits

        """
        if not isinstance(other, SequentialOrbits):
            return False
        # else:
        tests = [
            getattr(self, n) == getattr(other, n)
            for n in ("_ro", "_vo", "_zo", "_data")
        ]
        tests.append(all(self._solarmotion == other._solarmotion))
        # tests.append(self._pot == other._pot)
        if all(tests):
            return True
        return False

    # /def

    # +------------------------------+
    # Serialize:

    def __getstate__(self):
        """__getstate__."""
        return self.__dict__

    # /def

    def __setstate__(self, state):
        """__setstate__."""
        self.__dict__ = state

    # /def

    # def __getstate__(self):
    #     return (self._data,
    #             self._ro, self._vo, self._zo, self._solarmotion,
    #             self._pot)
    # # /def
    #
    # def __setstate__(self, state):
    #     print(state)
    #     _data, ro, vo, zo, solarmotion, pot = state
    #     self._data = _data
    #     self._ro = ro
    #     self._vo = vo
    #     self._zo = zo
    #     self._solarmotion = solarmotion
    #     self._pot = pot
    # # /def

    def dump(self, fname, protocol=None, *, fopt="b", fix_imports=True):
        """dump."""
        _dump(
            self, fname, protocol=protocol, fopt=fopt, fix_imports=fix_imports
        )

    # /def

    def save(self, fname, protocol=None, *, fopt="b", fix_imports=True):
        """save."""
        self.dump(fname, protocol=protocol, fopt=fopt, fix_imports=fix_imports)

    # /def

    @staticmethod
    def load(
        fname, *, fopt="b", fix_imports=True, encoding="ASCII", errors="strict"
    ):
        """load."""
        self = _load(
            fname,
            fopt=fopt,
            fix_imports=fix_imports,
            encoding=encoding,
            errors=errors,
        )
        return self

    # /def

    def copy(self):
        """copy."""
        instance = self.__class__(
            vxvv=self._data,
            pot=self._pot,
            solarmotion=self._solarmotion,
            ro=self._ro,
            vo=self._vo,
            zo=self._zo,
        )
        return instance

    # /def

    # +---------- time iteration ----------+
    def iterator(self, time=None):
        """iterator."""
        if time is None:
            timeiter = ((0, self._data.t0),)
        else:
            timeiter = self._data.iterator(time)
        return timeiter

    # /def

    # +------------------------------+
    # Add Orbits

    def _check_curr_orbit_integrated(self):
        """_check_curr_orbit_integrated.

        TODO check for None in _times instead

        """
        # Checking integrated previous orbit (has non-zero time domain)
        try:
            len(self.time)
        except TypeError:  # time is scalar, not a range
            warn(
                "Never integrated previous Orbit. Keeping current.",
                category=integrationWarning,
            )
            return False
        return True

    # /def

    def add_orbit(self, drct="previous", pot=None, vxvv=None, **kw):
        """Add orbit.

        Parameters
        ----------
        drct: str
            'previous': add orbit from previous orbit
            'forward': add orbit going forward in time
            'backward': add orbit going backward in time
        pot: Potential
        vxvv: list

        Returns
        -------
        self

        """
        if drct in ("prev", "previous"):
            self.add_orbit_from_prev(pot=pot, vxvv=vxvv, **kw)
        elif drct in ("for", "forward", "forwards"):
            self.add_forward_orbit(pot=pot, vxvv=vxvv, **kw)
        elif drct in ("back", "backward", "backwards"):
            self.add_backward_orbit(pot=pot, vxvv=vxvv, **kw)
        else:
            raise ValueError(drct)

        return self

    # /def

    def add_orbit_from_prev(self, pot=None, vxvv=None, **kw):
        """Add orbit from prev conditions.

        Parameters
        ----------
        pot: Potential
        vxvv: list

        Returns
        -------
        self

        """
        # Checking integrated previous orbit (has non-zero time domain)
        if not self._check_curr_orbit_integrated():
            return

        # Forward or Backwards?
        drct = self._data.direction
        if drct == "forward":
            oldt = self._bounds[1]  # getting current time

            if vxvv is None:  # use most recent orbit
                orbit = Orbit(vxvv=self.o.SkyCoord(oldt))
            elif isinstance(vxvv, Orbit):
                orbit = vxvv
            else:
                orbit = Orbit(vxvv=vxvv, **kw)

            self._data.append(oldt, orbit)

        else:
            oldt = self._bounds[0]  # getting current time

            if vxvv is None:  # use most recent orbit
                orbit = Orbit(vxvv=self.o.SkyCoord(oldt))
            elif isinstance(vxvv, Orbit):
                orbit = vxvv
            else:
                orbit = Orbit(vxvv=vxvv, **kw)

            self._data.prepend(oldt, orbit)

        # potential
        if pot is not None:
            self._pot = pot

        return self

    # /def

    def add_forward_orbit(self, pot=None, vxvv=None, **kw):
        """add_forward_orbit.

        Parameters
        ----------
        pot
            potential
        vxvv
            initial phase parameters
        kw
            into Orbit

        Returns
        -------
        self

        """
        # Checking integrated previous orbit (has non-zero time domain)
        if not self._check_curr_orbit_integrated():
            return

        ind = self._data._indmax  # forward orbit key
        oldt = self._data.tmax

        # adding orbit
        if vxvv is None:  # use forwardest orbit
            orbit = Orbit(vxvv=self._orbits[ind].SkyCoord(oldt))
        elif isinstance(vxvv, Orbit):
            orbit = vxvv
        else:
            orbit = Orbit(vxvv=vxvv, **kw)

        self._data.append(oldt, orbit)
        self._data.direction = "forward"

        # potential
        # need to check to switch potential b/c can be switching location
        if pot is not None:  # set new default potential
            self._pot = pot
        else:  # get potential from forwardest orbit
            self._pot = self._orbits[ind]._pot

        return self

    # /def

    def add_backward_orbit(self, pot=None, vxvv=None, **kw):
        """Start a backward orbit.

        Parameters
        ----------
        pot: Potential, optional
            change the potential
        vxvv: array, optional
            change initial conditions for this orbit

        Returns
        -------
        self

        """
        ind = self._data._indmin  # backward orbit key
        oldt = self._data.tmin

        # adding orbit
        if vxvv is None:  # use backwardest orbit
            orbit = Orbit(vxvv=self._orbits[ind].SkyCoord(t=oldt))
        elif isinstance(vxvv, Orbit):  # add orbit as is
            orbit = vxvv
        else:
            orbit = Orbit(vxvv=vxvv, **kw)  # standard add orbit

        self._data.prepend(oldt, orbit)
        self._data.direction = "backward"

        # potential
        # need to check to switch potential b/c can be switching location
        if pot is not None:  # set new default potential
            self._pot = pot
        else:  # get potential from backwardest orbit
            self._pot = self._orbits[ind]._pot

        return self

    # /def

    # +------------------------------+
    # Integrating

    def integrate(
        self,
        t,
        pot=None,
        method="symplec4_c",
        dt=None,
        fromprevt=True,
        _print=False,
    ):
        """Galpy Orbit integrate wrapper.

        potential always should be the same

        In Galpy, the starting time value does not matter

        sets the time based off integration time, so don't just do (0,...)
        actually set the correct time

        Parameters
        ----------
        fromprevt:
            adjust time so that integrating from the end of the last orbit

        Returns
        -------
        self

        TODO
        ----
        assume that anything without times is in Gyr

        """
        # Adjust time
        if not issubclass(t.__class__, u.Quantity):
            t = t * u.Gyr
        if self._data.direction == "backward":
            if t[1] > t[0]:
                t = t[::-1]
        if fromprevt:
            t = self._data.tref + (t - t[0])

        # Potential
        if pot is not None:  # keep potential current
            self._pot = pot  # store potential
        else:  # potential already set
            pot = self._pot  # set potential

        # Integrating
        if _print:
            print(f"integrating: {t[0]:.2f} : {t[-1]:.2f} : ({len(t)})")
        self.o.integrate(t, pot, method=method, dt=dt)
        if _print:
            print(f"integrated")  # TODO make write on previous line

        # Cleanup
        drct = self._data.direction
        if drct == "forward":
            self.time = t
            self._bounds = [t[0], t[-1]]
        else:
            self.time = t[::-1]
            self._bounds = [t[-1], t[0]]

        return self

    # /def

    # +---------- SKYCOORD ----------+
    def SkyCoord(
        self,
        t=None,
        ind=None,
        *args,
        frame: Optional[SkyCoord]=None,
        return_t: bool=False,
        T: bool=False,
        **kw,
    ):
        """SkyCoord.

        Parameters
        ----------
        t: array_like
            time
        ind: int
            the index of orbits within each segment to take skycoord

        .. note::
            args and kw do nothing right now

        Returns
        -------
        SkyCoord
        ts
            if return_t is True

        """
        # output frame
        if frame is None:  # inherit current frame
            frame = self.o.SkyCoord()

        # astropy bug. not the same galcen_v_sun
        if hasattr(frame, "galcen_v_sun"):
            frame_galcen_v_sun = frame.galcen_v_sun
        else:
            frame_galcen_v_sun = None  # will update later

        # taking care of index options, slices & lists are naturally supported
        if ind in (None, Ellipsis):
            ind = slice(None)

        scs, ts = [], []  # initializing SkyCoord, time arrays

        # iterating through orbits
        for i, time in self.iterator(t):  # orbit segment, time w/in segment
            if not isinstance(i, int):
                raise ValueError("i not int")

            try:  # getting at times
                sc = self._orbits[i][ind].SkyCoord(time).transform_to(frame)
            except ValueError:  # not integrated, get initial conditions
                sc = self._orbits[i][ind].SkyCoord().transform_to(frame)
            except IndexError:
                sc = self._orbits[i].SkyCoord(time).transform_to(frame)

            # astropy bug. not the same, even when should be
            if hasattr(sc, "galcen_v_sun"):
                if frame_galcen_v_sun is None:  # b/c no provided frame
                    frame_galcen_v_sun = sc.galcen_v_sun
                sc.galcen_v_sun = frame_galcen_v_sun

            scs.append(sc.T)
            ts.append(time)
        # /for

        try:
            scs = concatenate(scs)
        except ValueError:  # only 1 point
            scs = scs[0]

        if not T:
            scs = scs.T

        if return_t:
            return scs, astrarray(ts).flatten()
        else:
            return scs

        # if return_t:  # return time as well as SkyCoord
        #     if T is True:  # untranspose
        #         return concatenate(scs), astrarray(ts).flatten()
        #     else:
        #         return concatenate(scs).T, astrarray(ts).flatten()
        # else:  # just return SkyCoord
        #     if T is True:
        #         return concatenate(scs)
        #     else:
        #         return concatenate(scs).T

    # /def

    def getOrbit(self):
        """Get Orbit.

        concatenate along time axis

        """
        orbs = [self._orbits[i].getOrbit() for i in range(self.numsegs)]
        return np.concatenate(orbs, axis=1)

    # /def

    ###########################################################################
    # Plotting

    def _selfinstr(self, s):
        if not isinstance(s, str):
            raise TypeError(f"{s} is not <str>")
        if s.startswith("self."):
            return s
        else:
            return "self." + s

    # /def

    def _attrprep(self, x, frame=None):
        if isinstance(x, str):
            return eval(f"{self._selfinstr(x)}")
        else:
            return x

    # /def

    @mpl_decorator(funcdoc=plt.hist)
    def hist(self, attr, bins=10, **kw):
        """Histogram.

        Parameters
        ----------
        attr: str or array
        bins:
            the bins for plt.hist

        """
        if "label" in kw:
            kw["label"] = str(kw["label"])
        # print(self._attrprep(attr))
        try:
            res = plt.hist(self._attrprep(attr), bins=bins, **kw)
        except AttributeError:
            res = plt.hist(self._attrprep(attr + "._data"), bins=bins, **kw)
        return res

    # /def

    @mpl_decorator(funcdoc=plt.scatter)
    def plot_orbit(
        self,
        d1="ra",
        d2="dec",
        t="full",
        ind=None,
        label=None,
        frame=None,
        **kw,
    ):
        """Plot the orbits.

        # TODO support changing representation_type

        """
        label = None if label in (None, False) else f"orbit {self.orbind}"
        # origin = self.SkyCoord(t=0*u.s, frame=frame)
        orbit = self.SkyCoord(
            t=t, ind=ind, T=True, frame=frame
        )  # TODO not SkyCoord method?

        xlabel = kw.pop("xlabel", d1)
        ylabel = kw.pop("ylabel", d2)

        # plt.scatter(getattr(origin, d1), getattr(origin, d2), label=label,
        #             s=20, c='r')
        res = plt.plot(
            getattr(orbit, d1),
            getattr(orbit, d2),
            label=label,
            xlabel=xlabel,
            ylabel=ylabel,
            **kw,
        )
        return res

    # /def

    def plotOrbit(
        self,
        d1="ra",
        d2="dec",
        t="full",
        ind=None,
        label=None,
        frame=None,
        **kw,
    ):
        """Plot the orbits."""
        return self.plot_orbit(
            d1=d1, d2=d2, t=t, ind=ind, label=label, frame=frame, **kw
        )

    # /def

    @mpl_decorator()
    def plot(self, x, y, pltype="plot", **kw):
        """General plot function.

        Parameters
        ----------
        x: str
        y: str
        pltype: str
            pltype for plt.plot

        TODO
        ----
        support changing frame, representation_type

        """
        # x & y
        x = self._attrprep(x)
        y = self._attrprep(y)

        if pltype == "scatter":
            line = plt.scatter(x, y, **kw)
        elif pltype == "errorbar":
            xerr = self._attrprep(kw.pop("xerr", None))
            yerr = self._attrprep(kw.pop("yerr", None))

            line = plt.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)
        else:
            line = plt.plot(x, y, pltype=pltype, **kw)

        return line

    # /def

    @mpl_decorator(funcdoc=plt.scatter)
    def scatter(self, x, y, **kw):
        """Plot with pltype='scatter'."""
        return self.plot(x, y, pltype="scatter", **kw)

    # /def

    @mpl_decorator(funcdoc=plt.errorbar)
    def errorbar(self, x, y, xerr=None, yerr=None, **kw):
        """Plot errorbar.

        Parameters
        ----------
        x: str
        y: str
        xerr
        yerr
        key-word arguments

        """
        return self.plot(x, y, pltype="errorbar", xerr=xerr, yerr=yerr, **kw)

    # /def

    #############################################
    # Plot Helper Methods

    @classmethod
    def plotoptions(cls):
        """Print plotting functions & their signatures.

        Info
        ----
        General
            .plot(x, y, pltype, )
            .scatter(x, y, )
            .errorbar(x, y, xerr, yerr, )
            .hist(attr, bins, )

        Specific
            .plotOrbit(label, s):
                plot an orbit

        """
        print(cls.plothelp.__doc__)

    # /def

    @classmethod
    def printcmaps(cls, **kw):
        """Colormap options.

        Perceptually Uniform Sequential:
            'viridis', 'plasma', 'inferno', 'magma'
        Sequential:
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
        Sequential:
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper'
        Diverging:
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
        Qualitative:
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c'
        Miscellaneous:
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'

        """
        print(cls.printcmaps.__doc__, **kw)

    # /def
