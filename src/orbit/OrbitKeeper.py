# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : OrbitKeeper
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""OrbitKeeper."""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# GENERAL
import numpy as np


# PROJECT-SPECIFIC
from .exceptions import warn, integrationWarning
from ._indexdict import IndexDict


#############################################################################
# CODE
#############################################################################


class OrbitKeeper(object):
    """OrbitKeeper.

    store orbit and time information for SequentialOrbits

    indices are integers with corresponding orbit in _pworbits
    The location is the time organization

    """

    # def __new__(cls, *args):
    #     """OrbitKeeper __new__ method
    #     used for copying OrbitKeeper
    #     """
    #     print(args)
    #     if len(args) == 1:  # orbit first argument
    #         print('inargs')
    #         orbit = args[0]  # orbitkeeper
    #
    #         # for copying
    #         if isinstance(orbit, OrbitKeeper):
    #             # instance = cls.__new__(cls)
    #             # instance.__init__(orbit._orbits[0], orbit.t0)
    #             return orbit[:]  # using copying
    #         else:
    #             raise Exception('must be OrbitKeeper')
    #     else:
    #         print('outargs')
    #         self = super().__new__(cls)
    #         return self
    # # /def

    def __init__(self, orbit, t0):
        """Instantiate OrbitKeeper.

        sets:
            _ind: current orbit index
            _direction: 'forward'
            t0: start time for starting orbit
            tref: start time of current orbit
            _times: IndexDict of orbit evaluation times
            _bounds: IndexDict of orbit time bounds
            _orbits: IndexDict of orbit instances
        """
        self._ind = 0  # current orbit index
        self._direction = "forward"  # direction

        # The Time
        self.t0 = t0  # start time of 1st orbit
        self.tref = t0  # start time of current orbit

        self._times = IndexDict([(0, None)])  # set(0, None)
        self._bounds = IndexDict([(0, None)])  # set(0, None)

        # The Orbits
        if isinstance(orbit, IndexDict):
            self._orbits = orbit
        else:
            self._orbits = IndexDict([(0, orbit)])  # set(0, orbit)
            # self._orbits.direction = self.direction

        return

    # /def

    # @classmethod
    # def loadOrbitKeeper(orbits):  # TODO
    #     self._ind = orbits._ind
    #     self.direction
    #     self.t0 = orbits[orbits.minkey].time()[0]
    #     self._times =

    def __getitem__(self, ind):
        """Get a slice of OrbitKeeper.

        for slice objects, can only accept slice(None)
        TODO: view vs copy?

        """
        if isinstance(ind, slice) and not ind == slice(None):
            raise ValueError

        ok = OrbitKeeper(None, self.t0)
        ok._ind = self._ind
        ok._direction = self._direction
        ok.tref = self.tref
        ok._times = self._times[ind]
        ok._orbits = self._orbits[ind]

        # ensure never subindex _bounds
        if np.isscalar(ind) or ind in (Ellipsis, None, slice(None)):
            ok._bounds = self._bounds[ind]
        elif not isinstance(
            ind, tuple
        ):  # not multidimensional  # TODO improve
            ok._bounds = self._bounds[ind]
        else:
            ok._bounds = self._bounds[ind[0]]

        return ok

    def __eq__(self, other):
        """Test for complete equality."""
        tests = [
            (getattr(self, n) == getattr(other, n))
            for n in ("_ind", "_direction", "t0", "tref")
        ]
        tests.append(self._bounds == other._bounds)
        tests.append(self._orbits == other._orbits)
        tests.append(all(self._times == other._times))

        if all(tests):
            return True
        return False

    def copy(self):
        """Copy method."""
        ok = OrbitKeeper(None, self.t0)
        ok._ind = self._ind
        ok._direction = self._direction
        ok.tref = self.tref
        ok._times = self._times
        ok._orbits = self._orbits
        ok._bounds = self._bounds
        return ok

    # +---------- time ----------+
    @property
    def time(self):
        """Current time."""
        return self._times[self._ind]

    # /def

    @time.setter
    def time(self, value):
        """Set current time."""
        self._times[self._ind] = value

    # /def

    @property
    def alltimes(self):
        """All times."""
        return self._times

    # /def

    # +---------- bounds ----------+
    @property
    def bounds(self):
        """Current bounds."""
        return self._bounds[self._ind]

    # /def

    @bounds.setter
    def bounds(self, value):
        """Current bounds."""
        self._bounds[self._ind] = value

    # /def

    @property
    def allbounds(self):
        """Current bounds."""
        return self._bounds

    # /def

    # +---------- tmax/min ----------+
    @property
    def tmax(self):
        """Maximum time."""
        # last of time-ordered bounds, take forward bound
        return self._bounds.valueslist(-1)[-1]

    # /def

    @property
    def tmin(self):
        """Minimum time."""
        # last of time-ordered bounds, take forward bound
        return self._bounds.valueslist(0)[0]

    # /def

    # +---------- indmax/min ----------+
    @property
    def _indmax(self):
        """key/index at end."""
        # last of time-ordered bounds, take forward bound
        return self._bounds.keyslist(-1)

    # /def

    @property
    def _indmin(self):
        """key/index at start."""
        # last of time-ordered bounds, take forward bound
        return self._bounds.keyslist(0)

    # /def

    # +---------- orbit ----------+
    @property
    def orbit(self):
        """Current orbit."""
        return self._orbits[self._ind]

    # /def

    @property
    def allorbits(self):
        """All orbit segments."""
        return self._orbits

    # /def

    # +---------- direction ----------+

    @property
    def direction(self):
        """direction."""
        return self._direction

    # /def
    @direction.setter
    def direction(self, value):
        self._direction = value
        # self._orbits.direction = value

    # /def

    # +---------- keys / keyslist ----------+
    def keys(self):
        """keys."""
        return self._times.keys()

    # /def

    def keyslist(self):
        """keyslist."""
        return self._times.keyslist()

    # /def

    # +---------- remove times, bounds, orbits ----------+
    def remove(self, key):
        """remove."""
        print(f"removing orbit #{key}")
        del self._times[key]
        del self._bounds[key]
        del self._orbits[key]

    # /def

    # +---------- pre/append ----------+
    def prepend(self, tstart, orbit):
        """prepend."""
        if not (
            self._ind
            == self._times.maxkey
            == self._bounds.maxkey
            == self._orbits.maxkey
        ):
            raise IndexError("some index went wrong")
        self._ind += 1  # increment orbit index
        # setting
        self._times.prepend(None)
        self._bounds.prepend(None)
        self._orbits.prepend(orbit)
        # cleanup
        self.tref = tstart  # ref time -> current time

    # /def

    def append(self, tstart, orbit):
        """Append."""
        if not (
            self._ind
            == self._times.maxkey
            == self._bounds.maxkey
            == self._orbits.maxkey
        ):
            raise IndexError("some index went wrong")
        self._ind += 1  # increment orbit index
        # setting
        self._times.append(None)
        self._bounds.append(None)
        self._orbits.append(orbit)
        # cleanup
        self.tref = tstart  # ref time -> current time

    # /def

    # +---------- time iteration ----------+
    def iterator(self, time):
        """Make time iterator.

        Parameters
        ----------
        time: str, list, int
            the times at which to evaluate the orbits
            options:
            - str: 'full'
            - list: list of times
            - integer: the orbit index

        Returns
        -------
        t: list of lists
            the times

        """
        # time is integer index of an orbit
        if isinstance(time, (int, np.integer)):
            i = time  # time is index
            return ((i, self._times[i]),)

        # the time
        t = []

        if isinstance(time, str):
            if time not in ("full", "start", "end"):
                ValueError("f{time} not 'full', 'start', 'end'")
            # time shortcuts
            elif time == "full":
                tind = slice(None)  # time array at evaluation
            elif time == "start":
                tind = 0
            elif time == "ending":
                tind = -1

            # break up the time into the correct ranges for the orbit segments
            bnds = self._bounds.itemslist()

            for i, bnd in bnds:
                if bnd is None:  # bound is not integrated
                    warn("some segment is not integrated", integrationWarning)
                    continue  # skip to next bounds

                t.append((i, self._times[i][tind]))

            return t

        # break up the time into the correct ranges for the orbit segments
        bnds = self._bounds.itemslist()
        # bnds = [b for b in bnds if b is not None]  # skip not integrated

        # iterating over orbit segments
        # i is not in numerical order
        for i, bnd in bnds:
            if bnd is None:  # bound is not integrated
                warn("some segment is not integrated", integrationWarning)
                continue  # skip to next bounds

            ind = (bnd[0] <= time) & (time < bnd[1])  # w/in bounds

            # single bool (b/c any(True) raises an error)
            if isinstance(ind, (bool, np.bool_)):
                if ind:  # True
                    t.append((i, time))
            # bool array
            elif any(ind):  # w/in bounds ?
                t.append((i, time[ind]))

        # figuring out if evaluating outside of any bounds
        indbnd = np.array(
            [(b[0] <= time) & (time < b[1]) for i, b in bnds]
        ).sum(
            axis=0, dtype=bool
        )  # time in any bins
        splitat = np.where(np.diff(indbnd) != 0)[0] + 1  # split at bin edges

        hasvals = np.array(
            list(map(any, np.split(indbnd, splitat)))
        )  # in /out of bin
        indsplit = np.split(np.arange(len(time)), splitat)  # bin edges indices

        for ind, used in zip(indsplit, hasvals):  # iterating through bins
            if not used:  # time outside orbit integrations
                warn(
                    "Not including between {}:{}".format(
                        time[ind[0]], time[ind[-1]]
                    ),
                    integrationWarning,
                )

        return t

    # /def
