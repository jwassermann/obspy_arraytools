#!/usr/bin/env python
"""
Seismic array class.
:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from collections import defaultdict
import collections
import copy
import math
import tempfile
import os
import shutil
import warnings
import sys

import numpy as np
import scipy as sp
from obspy.core.util import SCIPY_VERSION
if SCIPY_VERSION < [1,14]:
    from scipy.integrate import cumtrapz
else:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects

from obspy.core import Trace, Stream, UTCDateTime
from obspy.core.event.event import Event
from obspy.core.event.origin import Origin
from obspy.core.inventory import Inventory
from obspy.core.util import AttribDict
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.imaging.cm import obspy_sequential
from obspy.signal.headers import clibsignal
from obspy.signal.invsim import cosine_taper
from obspy.signal.util import util_geo_km, next_pow_2
from obspy.taup import TauPyModel

from .beamforming_result import BeamformerResult
from .beamforming_result import plot_array_analysis
from .array_rotation_strain import \
    array_rotation_strain



def _get_stream_offsets(stream, stime, etime):
    """
    Calculates start and end offsets relative to stime and etime for each
    trace in stream in samples.

    :type stime: :class:`~obspy.core.utcdatetime.UTCDateTime`
    :param stime: Start time
    :type etime: :class:`~obspy.core.utcdatetime.UTCDateTime`
    :param etime: End time
    :returns: start and end sample offset arrays
    """
    stream = stream.sort()
    spoint = np.empty(len(stream), dtype=np.int32, order="C")
    epoint = np.empty(len(stream), dtype=np.int32, order="C")
    for i, tr in enumerate(stream):
        #if tr.stats.starttime > stime:
        if np.abs(tr.stats.starttime - stime) > tr.stats.sampling_rate:
            msg = "Specified stime %s is smaller than starttime %s " \
                  "in stream"
            raise ValueError(msg % (stime, tr.stats.starttime))
       # if tr.stats.endtime < etime:
        if np.abs(tr.stats.endtime - etime) > tr.stats.sampling_rate:
            msg = "Specified etime %s is bigger than endtime %s in stream"
            raise ValueError(msg % (etime, tr.stats.endtime))
        # now we have to adjust to the beginning of real start time
        spoint[i] = int(
            (stime - tr.stats.starttime) * tr.stats.sampling_rate + .5)
        epoint[i] = int(
            (tr.stats.endtime - etime) * tr.stats.sampling_rate + .5)
    return spoint, epoint


class SeismicArray(object):
    """
    Class representing a seismic (or other) array.

    The SeismicArray class is a named container for an
    :class:`~obspy.core.inventory.Inventory` containing the components
    making up the array along with methods for array processing. It does not
    contain any seismic data. The locations of the array components
    (stations or channels) must be set in the respective objects (for an
    overview of the inventory system, see :mod:`~obspy.core.inventory`).
    While the inventory must be composed of
    :class:`~obspy.core.inventory.channel.Channel` objects,
    they do not need to represent actual seismic stations; their only required
    attribute is the location information.

    :param name: Array name.
    :type name: str
    :param inventory: Inventory of stations making up the array.
    :type inventory: :class:`~obspy.core.inventory.Inventory`

    .. rubric:: Basic Usage

    >>> from obspy.core.inventory import read_inventory
    >>> from obspy.signal.array_analysis import SeismicArray
    >>> inv = read_inventory('http://examples.obspy.org/agfainventory.xml')
    >>> array = SeismicArray('AGFA', inv)
    >>> print(array)
    Seismic Array 'AGFA' with 5 Stations, aperture: 0.06 km.

    .. rubric:: Coordinate conventions:

    * Right handed
    * X positive to east
    * Y positive to north
    * Z positive up
    """

    def __init__(self, name, inventory):
        """
        :param name: Name of this Array.
        :type name: str
        :param inventory: Inventory containing the stations used for the array.
        :type inventory: :class:`~obspy.core.inventory.inventory.Inventory`
        """
        if not isinstance(name, str):
            raise TypeError("Name must be a string.")
        self.name = name
        if not isinstance(inventory, Inventory):
            raise TypeError("Inventory must be an ObsPy Inventory.")
        self.inventory = copy.deepcopy(inventory)

    def __str__(self):
        """
        Pretty representation of the array.
        """
        if self.inventory is None:
            return "Empty seismic array '{}'".format(self.name)
        ret_str = "Seismic Array '{name}' with ".format(name=self.name)
        ret_str += "{count} Stations, ".format(count=len(self.geometry))
        ret_str += "Aperture: {aperture:.2f} km.".format(
            aperture=self.aperture)
        return ret_str

    def inventory_cull(self, stream):
        """
        Shrink array inventory to channels present in the given stream.

        Permanently remove from the array inventory all entries for stations or
        channels that do not have traces in the given
        :class:`~obspy.core.stream.Stream` st. This may be useful e.g. if data
        is not consistently available from every channel in an array, or if a
        web service has returned many more inventory objects than required. The
        method selects channels based on matching network, station, location
        and channel codes to the ones given in the trace headers. Furthermore,
        if a time range is specified for a channel, it will only be kept if it
        matches the time span of its corresponding trace.

        If you wish to keep the original inventory, make a copy first:

        >>> from copy import deepcopy #doctest: +SKIP
        >>> original_inventory = deepcopy(array.inventory) #doctest: +SKIP

        :param stream: :class:`~obspy.core.stream.Stream` to which the array
            inventory should correspond.
        :type stream: :class:`~obspy.core.stream.Stream`
        """
        inv = self.inventory
        stream = stream.sort()
        # check what station/channel IDs are in the data
        stations_present = list(set(tr.id for tr in stream))
        # delete all channels that are not represented
        for k, netw in reversed(list(enumerate(inv.networks))):
            for j, stn in reversed(list(enumerate(netw.stations))):
                for i, cha in reversed(list(enumerate(stn.channels))):
                    if ("{}.{}.{}.{}".format(netw.code, stn.code,
                                             cha.location_code, cha.code)
                            not in stations_present):
                        del stn.channels[i]
                        continue
                        # Also remove if it doesn't cover the time of the
                        # trace:
                    for tr in stream.select(network=netw.code,
                                            station=stn.code,
                                            location=cha.location_code,
                                            channel=cha.code):
                        if not cha.is_active(starttime=tr.stats.starttime,
                                             endtime=tr.stats.endtime):
                            del stn.channels[i]
                stn.total_number_of_channels = len(stn.channels)
                # no point keeping stations with all channels removed
                if len(stn.channels) == 0:
                    del netw.stations[j]
            # no point keeping networks with no stations in them:
            if len(netw.stations) == 0:
                del inv.networks[k]
        # check total number of channels now:
        contents = inv.get_contents()
        if len(contents['channels']) < len(stations_present):
            # Inventory is altered anyway in this case.
            warnings.warn('Inventory does not contain information for all '
                          'traces in stream.')
        self.inventory = inv

    def plot(self, projection="local", show=True, **kwargs):
        """
        Plot the geographical layout of the array.

        Shows all the array's stations as well as it's geometric center and its
        center of gravity.
        >>> from obspy.core.inventory import read_inventory
        >>> from obspy.signal.array_analysis import SeismicArray
        >>> inv = read_inventory('http://examples.obspy.org/agfainventory.xml')
        >>> array = SeismicArray('AGFA', inv)
        >>> array.plot()

        .. plot::

            from obspy.core.inventory import read_inventory
            from obspy.signal.array_analysis import SeismicArray
            inv = read_inventory('http://examples.obspy.org/agfainventory.xml')
            array = SeismicArray('AGFA', inv)
            array.plot()

        :type projection: str, optional
        :param projection: The map projection. Currently supported are:

            * ``"global"`` (Will plot the whole world.)
            * ``"ortho"`` (Will center around the mean lat/long.)
            * ``"local"`` (Will plot around local events)

            Defaults to ``"local"``
        :type show: bool
        :param show: Whether to show the figure after plotting or not. Can be
            used to do further customization of the plot before showing it.

        All other keyword arguments are passed to the
        :meth:`obspy.core.inventory.inventory.Inventory.plot` method.
        """
        # Piggy-back on the inventory plotting. Currently requires basemap.
        fig = self.inventory.plot(projection=projection, show=False, **kwargs)
                                  #method="basemap", **kwargs)
        bmap = fig.bmap

        path_effects = [patheffects.withStroke(linewidth=3,
                                               foreground="white")]

        grav = self.center_of_gravity
        x, y = bmap(grav["longitude"], grav["latitude"])
        #x = grav["longitude"]; y=grav["latitude"]
        bmp = fig.axes[0]
        bmp.scatter(x, y, marker="x", c="blue", s=100, zorder=201,
                     linewidths=2)
        bmp.text(x, y, " Center of Gravity", color="blue", ha="left",
                     weight="heavy", zorder=200,
                     path_effects=path_effects)

        geo = self.geometrical_center
        x, y = bmap(geo["longitude"], geo["latitude"])
        #x = geo["longitude"]; y= geo["latitude"]
        bmp.scatter(x, y, marker="x", c="green", s=100, zorder=201,
                     linewidths=2)
        bmp.text(x, y, "Geometrical Center ", color="green", ha="right",
                     fontweight=900, zorder=200, path_effects=path_effects)

        bmp.set_title(str(self).splitlines()[0].strip())

        if show:
            plt.show()
        return fig

    def _get_geometry(self):
        """
        Return a dictionary of latitude, longitude and absolute height
        [km] for each component in the array inventory.

        For every component in the array inventory (channels if available,
        stations otherwise), a SEED ID string with the format
        'network.station.location.channel', leaving any unknown parts blank, is
        assembled. This is one key for the returned dictionary, while the value
        is a dictionary of the component's coordinates.

        :return A dictionary with keys: SEED IDs and values: dictionaries of
            'latitude', 'longitude' and 'absolute_height_in_km'.
        """
        if not self.inventory:
            return {}
        geo = {}

        # Using core.inventory.inventory.Inventory.get_coordinates() is not
        # really satisfactory: It doesn't return coordinates for inventories
        # that have stations but no channels defined.
        # Might be the case e.g. if using the array class for inventory of
        # sources.
        for network in self.inventory:
            for station in network:
                if len(station.channels) == 0:
                    # Using the full Seed ID string allows retrieving
                    # coordinates with geometry[trace.id] for other methods.
                    item_code = "{n}.{s}..".format(n=network.code,
                                                   s=station.code)
                    this_coordinates = \
                        {"latitude": float(station.latitude),
                         "longitude": float(station.longitude),
                         "absolute_height_in_km":
                         float(station.elevation) / 1000.0}
                    geo[item_code] = this_coordinates
                else:
                    for channel in station:
                        item_code = "{}.{}.{}.{}".format(network.code,
                                                         station.code,
                                                         channel.location_code,
                                                         channel.code)
                        this_coordinates = \
                            {"latitude": float(channel.latitude),
                             "longitude": float(channel.longitude),
                             "absolute_height_in_km":
                             float(channel.elevation - channel.depth) / 1000.0}
                        geo[item_code] = this_coordinates
        return geo

    @property
    def geometry(self):
        """
        A dictionary of latitude, longitude and absolute height [km] values
        for each item in the array inventory.
        For every component in the array inventory (channels if available,
        stations otherwise), a SEED ID string with the format
        'network.station.location.channel', leaving any unknown parts blank, is
        assembled. This is one key for the returned dictionary, while the value
        is a dictionary of the component's coordinates.
        :return A dictionary with keys: SEED IDs and values: dictionaries of
            'latitude', 'longitude' and 'absolute_height_in_km'.
        """
        return self._get_geometry()

    @property
    def geometrical_center(self):
        """
        Return the geometrical centre as dictionary.
        The geometrical centre is the mid-point of the maximum array extent in
        each direction.

        :return A dictionary with keys: latitude, longitude and
        absolute_height_in_km:
        """
        extent = self.extent
        return {
            "latitude": (extent["max_latitude"] +
                         extent["min_latitude"]) / 2.0,
            "longitude": (extent["max_longitude"] +
                          extent["min_longitude"]) / 2.0,
            "absolute_height_in_km":
            (extent["min_absolute_height_in_km"] +
             extent["max_absolute_height_in_km"]) / 2.0
        }

    @property
    def center_of_gravity(self):
        """
        Return the centre of gravity as a dictionary.
        The centre of gravity is calculated as the mean of the array stations'
        locations in each direction.

        :return A dictionary with keys: latitude, longitude and
        absolute_height_in_km:
        """
        lats, lngs, hgts = self._coordinate_values()
        return {
            "latitude": np.mean(lats),
            "longitude": np.mean(lngs),
            "absolute_height_in_km": np.mean(hgts)}

    @property
    def aperture(self):
        """
        Return the array aperture in kilometers.
        The array aperture is the maximum distance between any two stations in
        the array.

        :return Array aperture in kilometers.
        """
        distances = []
        geo = self.geometry
        for location, coordinates in geo.items():
            for other_location, other_coordinates in list(geo.items()):
                if location == other_location:
                    continue
                distances.append(gps2dist_azimuth(
                    coordinates["latitude"], coordinates["longitude"],
                    other_coordinates["latitude"],
                    other_coordinates["longitude"])[0] / 1000.0)
        return max(distances)

    @property
    def extent(self):
        """
        Dictionary of the array's minimum and maximum lat/long and elevation
        values.

        :return A dictionary with keys: min_latitude, max_latitude,
        min_longitude, max_longitude, min_absolute_height_in_km and
        max_absolute_height_in_km:
        """
        lats, lngs, hgt = self._coordinate_values()

        return {
            "min_latitude": min(lats),
            "max_latitude": max(lats),
            "min_longitude": min(lngs),
            "max_longitude": max(lngs),
            "min_absolute_height_in_km": min(hgt),
            "max_absolute_height_in_km": max(hgt)}

    def _coordinate_values(self):
        """
        Return the array geometry as simple lists of latitude, longitude and
        elevation.
        """
        geo = self.geometry
        lats, lngs, hgt = [], [], []
        for coordinates in list(geo.values()):
            lats.append(coordinates["latitude"]),
            lngs.append(coordinates["longitude"]),
            hgt.append(coordinates["absolute_height_in_km"])
        return lats, lngs, hgt

    def _get_geometry_xyz(self, latitude, longitude, absolute_height_in_km,
                          correct_3dplane=False):
        """
        Return the array geometry as each components's offset relative to a
        given reference point, in km.

        The returned geometry is a nested dictionary with each component's SEED
        ID as key and a dictionary of its coordinates as value, similar to that
        given by :attr:`~obspy.signal.array_analysis.SeismicArray.geometry`,
        but with different coordinate values and keys.

        To obtain the x-y-z geometry in relation to, for example, the center of
        gravity, use:

        >>> array = SeismicArray('', inv) # doctest: +SKIP
        >>> array._get_geometry_xyz(**array.center_of_gravity) # doctest: +SKIP

        :param latitude: Latitude of reference origin.
        :param longitude: Longitude of reference origin.
        :param absolute_height_in_km: Elevation of reference origin.
        :type correct_3dplane: bool
        :param correct_3dplane: Correct the returned geometry by a
            best-fitting 3D plane.
            This might be important if the array is located on an inclined
            slope.
        :return: The geometry of the components as dictionary, with coordinate
            keys of 'x', 'y' and 'z'.
        """
        geometry = {}
        for key, value in list(self.geometry.items()):
            x, y = util_geo_km(longitude, latitude, value["longitude"],
                               value["latitude"])
            geometry[key] = {
                "x": x,
                "y": y,
                "z": value["absolute_height_in_km"] - absolute_height_in_km
            }
        if correct_3dplane:
            geometry = self._correct_with_3dplane(geometry)
        return geometry

    def _get_timeshift_baz(self, sll, slm, sls, baz, latitude, longitude,
                           absolute_height_in_km, static3d=False, vel_cor=4.0, sec_km=False):
        """
        Returns timeshift table for the geometry of the current array, in
        kilometres relative to a given centre (uses geometric centre if not
        specified), and a pre-defined backazimuth. Returns nested dict of
        timeshifts at each slowness between sll and slm, with sls increment.

        :param sll: slowness x min (lower)
        :param slm: slowness x max (lower)
        :param sls: slowness step
        :param baz:  backazimuth applied
        :param latitude: latitude of reference origin
        :param longitude: longitude of reference origin
        :param absolute_height_in_km: elevation of reference origin, in km
        :param vel_cor: Correction velocity (upper layer) in km/s. May be given
            at each station as a dictionary with the station/channel IDs as
            keys (same as in self.geometry).
        :type static3d: bool
        :param static3d: a correction of the station height is applied using
            vel_cor the correction is done according to the formula:
            t = rxy*s - rz*cos(inc)/vel_cor
            where inc is defined by inc = asin(vel_cor*slow)
        :param sec_km: switch for slowness in and output. If true s/km values are expected
                       if false (default) than s/degree are used
        :return Dictionary with time differences relative to reference point.
            The station names are given in the keys and a lost of time
            differences for each slowness step is given in the items.
            The none key contains the absolute slowness value for each step.
        """
        geom = self._get_geometry_xyz(latitude, longitude,
                                      absolute_height_in_km)
        baz = math.pi * baz / 180.0
        if sec_km:
            KM_PER_DEG = 1.
        else:
            KM_PER_DEG = 111.1949


        time_shift_tbl = {}
        slownesses = np.arange(sll, slm + sls, sls)
        # we have to bring slowness (ray based coordinate system) and
        # baz (oposite direction) together -slowness
        time_shift_tbl[None] = slownesses

        for key, value in list(geom.items()):
            time_shifts = -slownesses / KM_PER_DEG * \
                    (value["x"] * math.sin(baz) + value["y"] * math.cos(baz))
            if static3d:
                try:
                    inc = np.arcsin(vel_cor * slownesses / KM_PER_DEG)
                except ValueError:
                    # if vel_cor given as dict:
                    inc = np.pi / 2.0
                try:
                    v = vel_cor[key]
                except TypeError:
                    # if vel_cor is a constant:
                    v = vel_cor
                time_shifts += value["z"] * np.cos(inc) / v
            time_shift_tbl[key] = time_shifts

        return time_shift_tbl

    def _get_timeshift(self, sllx, slly, sls, grdpts_x, grdpts_y,
                       latitude=None, longitude=None, absolute_height=None,
                       vel_cor=4., static3d=False, sec_km = False):
        """
        Returns timeshift table for the geometry of the current array, in
        kilometres relative to a given centre (uses geometric centre if not
        specified).

        :param sllx: slowness x min (lower)
        :param slly: slowness y min (lower)
        :param sls: slowness step
        :param grdpts_x: number of grid points in x direction
        :param grdpts_y: number of grid points in y direction
        :param latitude: latitude of reference origin
        :param longitude: longitude of reference origin
        :param absolute_height: elevation of reference origin, in km
        :param vel_cor: correction velocity (upper layer) in km/s
        :type static3d: bool
        :param static3d: a correction of the station height is applied using
            vel_cor the correction is done according to the formula:
            t = rxy*s - rz*cos(inc)/vel_cor
            where inc is defined by inc = asin(vel_cor*slow)
        :type sec_km: bool
        :param sec_km: switch for slowness in and output. If true s/km values are expected
                       if false (default) than s/degree are used
        :return 2D timeshift table for each station in the array. Each table
                gives the timeshift for all slowness_x and slowness_y
                combinations.
        """
        if sec_km == True:
            KM_PER_DEG = 1.
        else:
            KM_PER_DEG = 111.1949

        if any([_i is None for _i in [latitude, longitude, absolute_height]]):
            #latitude = self.geometrical_center["latitude"]
            #longitude = self.geometrical_center["longitude"]
            #absolute_height = self.geometrical_center["absolute_height_in_km"]
            latitude = self.center_of_gravity["latitude"]
            longitude = self.center_of_gravity["longitude"]
            absolute_height = self.center_of_gravity["absolute_height_in_km"]
        geom = self._get_geometry_xyz(latitude, longitude,
                                      absolute_height)

        geometry = self._geometry_dict_to_array(geom)
        if static3d:
            nstat = len(geometry)
            time_shift_tbl = np.empty((nstat, grdpts_x, grdpts_y),
                                      dtype="float32")
            for i in range(grdpts_x):
                sx = (sllx + i * sls) / KM_PER_DEG
                for j in range(grdpts_y):
                    sy = (slly + j * sls) / KM_PER_DEG
                    slow = np.sqrt(sx * sx + sy * sy)
                    if vel_cor * slow <= 1.:
                        inc = np.arcsin(vel_cor * slow)
                    else:
                        warnings.warn(
                            "Correction velocity smaller than apparent"
                            " velocity")
                        inc = np.pi / 2.
                    time_shift_tbl[:, i, j] = sx * geometry[:, 0] + sy * \
                        geometry[:, 1] + geometry[:, 2] * np.cos(inc) / vel_cor
            return time_shift_tbl
        # optimized version
        else:
            mx = np.outer(geometry[:, 0],
                          (sllx + np.arange(grdpts_x) * sls) / KM_PER_DEG)
            my = np.outer(geometry[:, 1],
                          (slly + np.arange(grdpts_y) * sls) / KM_PER_DEG)
            return np.require(
                mx[:, :, np.newaxis].repeat(grdpts_y, axis=2) +
                my[:, np.newaxis, :].repeat(grdpts_x, axis=1),
                dtype='float32')

    def vespagram(self, stream, event_or_baz, sll, slm, sls, starttime,
                  endtime, reference='center_of_gravity', method="DLS",
                  nthroot=1, static3d=False, vel_cor=4.0, wiggle_scale=1.0,
                  align=False, align_phase='P',
                  density_cmap=obspy_sequential, plot="wiggle", show=True):
        """
        :type stream: :class:`~obspy.core.stream.Stream`.
        :param stream: Stream containing all traces of the array.
        :type event_or_baz: float or :class:`~obspy.core.event.event.Event` or
            :class:`~obspy.core.event.origin.Origin`
        :param event_or_baz: Backazimuth for vespagram or event/origin object
            to calculate theoretical backazimuth from.
        :type sll: float
        :param sll: Lower slowness boundary.
        :type slm: float
        :param slm: Maximum slowness boundary.
        :type sls: float
        :param sls: Slowness step.
        :type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param starttime: Beginn of analysed sequence.
        :type endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param endtime: End of analysed sequence.
        :type reference: str or dict
        :param reference: Determines what is used as reference origin. Either
            ``'center_of_gravity'``, ``'geometrical_center'`` or a dictionary
            with keys ``'latitude'``, ``'longitude'``, ``'elevation'``
            (elevation in meters).
        :type method: str
        :param method: Method used to stack traces. Either 'DLS' or 'PWS':
                        DLS: delay and sum
                        PWS: phase-weighted stack
        :type nthroot: float
        :param nthroot: Root used in th summation process.
        :type static3d: bool
        :param static3d: a correction of the station height is applied using
            vel_cor the correction is done according to the formula:
            t = rxy*s - rz*cos(inc)/vel_cor
            where inc is defined by inc = asin(vel_cor*slow)
        :param vel_cor: float or dict
        :param vel_cor: correction velocity (upper layer) in km/s. May be given
            at each station as a dictionary with the station/channel IDs as
            keys (same as in self.geometry).
        :type wiggle_scale: float
        :param wiggle_scale: Relative scaling for wiggle plot
                            (unused for density plot).
        :type align: bool
        :param align: Whether traces should be aligned with theoretical arrival
                       given by phase_name.
        :type align_phase: str
        :param align_phase: Phase by which the traces are aligned. Mostly
                           conventional naming, for more information see:
                            docs.obspy.org/packages/obspy.taup.html
        :type density_cmap: :class:`~matplotlib.colors.Colormap`
        :param density_cmap: Colormap used for density plot.
        :type plot: str or None
        :param plot: Whether to create a plot or not. Can be either
            ``'wiggle'``, ``'density'``, or ``None``.
        :type show: bool
        :param show: Whether to open the plot (if any) interactively or not.
        :returns slow: Slowness with highest coherency.
                 beams: :class:`numpy.ndarray` with beam for each slowness step.
                 beam_max: Inde of beam with maximum coherency.
                 max_beam: :class:`numpy.ndarray` of beam with maximum coherency.
                 fig: :class:`matplotlib.pyplot.figure` of vespagramm.
        """
        stream = stream.sort()
        if reference == 'center_of_gravity':
            center_ = self.center_of_gravity
        elif reference == 'geometrical_center':
            center_ = self.geometrical_center
        elif isinstance(reference, dict):
            center_ = reference
            center_['absolute_height_in_km'] = (
                center_.pop('elevation') / 1000.0)
        else:
            msg = "Unrecognized value for 'reference' option: {}"
            raise ValueError(msg.format(reference))

        if isinstance(event_or_baz, Event):
            origin_ = event_or_baz.origins[0]
            baz = gps2dist_azimuth(
                center_['latitude'], center_['longitude'],
                origin_['latitude'], origin_['longitude'])[1]

        elif isinstance(event_or_baz, Origin):
            origin_ = event_or_baz
            baz = gps2dist_azimuth(
                center_['latitude'], center_['longitude'],
                origin_['latitude'], origin_['longitude'])[1]
        else:
            baz = float(event_or_baz)
            if align:
                msg = "For the align option an event has to be specified"
                raise ValueError(msg)

        if align:
            stream = self.align_phases(stream, origin_, align_phase)

        time_shift_table = self._get_timeshift_baz(
            sll, slm, sls, baz, latitude=center_['latitude'],
            longitude=center_['longitude'],
            absolute_height_in_km=center_['absolute_height_in_km'],
            static3d=static3d, vel_cor=vel_cor)

        slownesses = time_shift_table[None]

        slow, beams, beam_max, max_beam = self._vespagram_baz(
            stream, time_shift_table, starttime=starttime, endtime=endtime,
            method=method, nthroot=nthroot)

        if plot:
            if plot not in ('wiggle', 'density'):
                msg = "Unknown plotting option: '{!s}'".format(plot)
                raise ValueError(msg)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # XXX need to check that all sampling rates are equal!?
            sampling_rate = stream[0].stats.sampling_rate
            delta = 1 / sampling_rate
            npts = len(beams[0])
            t = np.linspace(0, npts*delta, npts)
            max_amp = np.max(np.abs(beams))
            scale = sls / max_amp
            scale *= wiggle_scale

            if plot == 'wiggle':
                for i_, beam in enumerate(beams):
                    if i_ == beam_max:
                        ax.plot(t, slownesses[i_] + scale * beams[i_], 'r',
                                zorder=1)
                    else:
                        ax.plot(t, slownesses[i_] + scale * beams[i_], 'k',
                                zorder=-1)
                ax.set_xlabel('Time [s]')
                ax.set_ylabel('slowness [s/deg]')
                ax.set_xlim(t[0], t[-1])
                data_minmax = ax.yaxis.get_data_interval()
                minmax = [min(slownesses[0], data_minmax[0]),
                          max(slownesses[-1], data_minmax[1])]
                ax.set_ylim(*minmax)

            elif plot == 'density':
                extent = (t[0] - delta * 0.5, t[-1] + delta * 0.5,
                          slownesses[0] - sls * 0.5,
                          slownesses[-1] + sls * 0.5)

                col = ax.imshow(np.flipud(beams), cmap=density_cmap,
                                interpolation="nearest", extent=extent,
                                aspect='auto')
                plt.colorbar(col)

            ax.set_ylabel('slowness [s/deg]')
            ax.set_xlabel('Time [s]')
            if show:
                plt.show()
        else:
            fig = None

        return slow, beams, beam_max, max_beam, fig

    def derive_rotation_from_array(self, stream, vp, vs, sigmau, latitude,
                                   longitude, absolute_height_in_km=0.0):
        """
        Wrapper for the Function
        :class:`~obspy.signal.array_analysis.array_rotation_strain.array_rotation_strain`.
        Returns rotations and strains calculations from array measurements in a
        structured way and pre-processes the input.
        :param stream: Stream containing all traces of the array.
        :type stream: :class:`~obspy.core.stream.Stream`.
        :param vp: P wave speed in the soil under the array (km/s).
        :param vs: S wave speed in the soil under the array Note - vp and vs
        may be any unit (e.g. miles/week), and this unit need not be related to
        the units of the station coordinates or ground motions, but the
        units of vp and vs must be the SAME because only their ratio is used.
        :param sigmau: Standard deviation (NOT VARIANCE) of ground noise,
            corresponds to sigma-sub-u in S95 lines above eqn (A5).
            NOTE: This may be entered as a scalar, vector, or matrix!
            * If sigmau is a scalar, it will be used for all components of all
              stations.
            * If sigmau is a 1D array of length Na, sigmau[i] will be the noise
              assigned to all components of the station corresponding to
              array_coords[i,:]
            * If sigmau is a 2D array of dimension  Na x 3, then sigmau[i,j] is
              used as the noise of station i, component j.
            In all cases, this routine assumes that the noise covariance
            between different stations and/or components is zero.
        :type sigmau: float or :class:`numpy.ndarray`
        :param latitude: Latitude of reference point.
        :param longitude: Longitude of reference point.
        :param absolute_height_in_km: Absolute height of reference point (km).
        :return: Rotated :class:`~obspy.core.stream.Stream`. and output of
         :class:`~obspy.signal.array_analysis.array_rotation_strain.array_rotation_strain`.
         which contains the rotation and strain parameters in a dictionary.
        """
        geo = self.geometry

        components = collections.defaultdict(list)
        for tr in stream:
            components[tr.stats.channel[-1].upper()].append(tr)

        # Sanity checks.
        if sorted(components.keys()) != ["E", "N", "Z"]:
            raise ValueError("Three components necessary.")

        for value in list(components.values()):
            value.sort(key=lambda x: "%s.%s" % (x.stats.network,
                                                x.stats.station))

        ids = [tuple([_i.id[:-1] for _i in traces]) for traces in
               list(components.values())]
        if len(set(ids)) != 1:
            raise ValueError("All stations need to have three components.")

        stats = [[(_i.stats.starttime.timestamp, _i.stats.npts,
                   _i.stats.sampling_rate)
                  for _i in traces] for traces in list(components.values())]
        s = []
        for st in stats:
            s.extend(st)

        if len(set(s)) != 1:
            raise ValueError("starttime, npts, and sampling rate must be "
                             "identical for all traces.")

        stations = ["%s.%s.%s.%s" % (_i.stats.network, _i.stats.station,
                                     _i.stats.location, _i.stats.channel)
                    for _i in list(components.values())[0]]
        for station in stations:
            if station not in geo:
                raise ValueError("No coordinates known for station '%s'" %
                                 station)
        # every third entry because every coordinate is represented three
        # times due to the three directions
        array_coords = self._geometry_dict_to_array(
            self._get_geometry_xyz(latitude,
                                   longitude,
                                   absolute_height_in_km))[::3]

        subarray = np.arange(len(geo)/3)
        # integer is needed for fancy indexing
        subarray = np.array(subarray, dtype=int)

        tr = []
        for _i, component in enumerate(["Z", "N", "E"]):
            comp = components[component]
            tr.append(np.empty((len(comp[0]), len(comp))))
            for _j, trace in enumerate(comp):
                tr[_i][:, _j][:] = np.require(trace.data, np.float64)

        rot = array_rotation_strain(subarray, tr[0], tr[1], tr[2], vp=vp,
                                    vs=vs, array_coords=array_coords,
                                    sigmau=sigmau)

        d1 = rot.pop("ts_w1")
        d2 = rot.pop("ts_w2")
        d3 = rot.pop("ts_w3")

        header = {"network": "XX", "station": "YY", "location": "99",
                  "starttime": list(components.values())[0][0].stats.starttime,
                  "sampling_rate":
                  list(components.values())[0][0].stats.sampling_rate,
                  "channel": "ROZ",
                  "npts": len(d1)}

        tr1 = Trace(data=d1, header=copy.copy(header))
        header["channel"] = "RON"
        header["npts"] = len(d2)
        tr2 = Trace(data=d2, header=copy.copy(header))
        header["channel"] = "ROE"
        header["npts"] = len(d3)
        tr3 = Trace(data=d3, header=copy.copy(header))

        return Stream(traces=[tr1, tr2, tr3]), rot

    def slowness_whitened_power(self, stream, frqlow, frqhigh,
                                prefilter=True, plots=(),
                                static3d=False,correct_3dplane=False, sec_km=False,array_response=False,verbose=False,
                                vel_corr=4.8, wlen=-1, wfrac=1.,
                                slx=(-10, 10), sly=(-10, 10), sls=0.5):
        """
        Slowness whitened power analysis.

        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :type static3d: bool
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :param array_response: Specify if array response should be used.
        :type array_response: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param wfrac: Fraction of the window not overlapping with other
                      windows.
        :type wfrac: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting options are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        return self._array_analysis_helper(stream=stream, method="SWP",
                                           frqlow=frqlow, frqhigh=frqhigh,
                                           prefilter=prefilter, plots=plots,
                                           static3d=static3d,correct_3dplane = correct_3dplane,sec_km=sec_km,
                                           array_r=array_response,verbose=verbose,
                                           vel_corr=vel_corr, wlen=wlen,
                                           wfrac=wfrac,
                                           slx=slx, sly=sly, sls=sls)

    def phase_weighted_stack(self, stream, frqlow, frqhigh,
                             prefilter=True, plots=(),verbose=False,
                             static3d=False,correct_3dplane=False,sec_km=False, array_response=False,
                             vel_corr=4.8, wlen=-1, wfrac=1., slx=(-10, 10),
                             sly=(-10, 10), sls=0.5):
        """
        Phase weighted stack analysis.

        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :type static3d: bool
        :param array_response: Specify if array response should be used.
        :type array_response: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type wfrac: float
        :param wfrac: Fraction of the window not overlapping with other
                      windows.
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting otions are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        return self._array_analysis_helper(stream=stream, method="PWS",
                                           frqlow=frqlow, frqhigh=frqhigh,
                                           prefilter=prefilter, plots=plots,
                                           static3d=static3d,correct_3dplane = correct_3dplane,sec_km=sec_km,
                                           array_r=array_response,verbose=verbose,
                                           vel_corr=vel_corr, wlen=wlen,
                                           wfrac=wfrac,
                                           slx=slx, sly=sly, sls=sls)

    def delay_and_sum(self, stream, frqlow, frqhigh,
                      prefilter=True, plots=(), static3d=False,correct_3dplane=False,sec_km=False,
                      array_response=False,verbose=False,
                      vel_corr=4.8, wlen=-1, wfrac=1., slx=(-10, 10),
                      sly=(-10, 10), sls=0.5):
        """
        Delay and sum analysis.
        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :type static3d: bool
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :param array_response: Specify if array response should be used.
        :type array_response: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param wfrac: Fraction of the window not overlapping with other
                      windows.
        :type wfrac: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting otions are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        return self._array_analysis_helper(stream=stream, method="DLS",
                                           frqlow=frqlow, frqhigh=frqhigh,
                                           prefilter=prefilter, plots=plots,
                                           static3d=static3d,correct_3dplane = correct_3dplane,sec_km=sec_km,
                                           array_r=array_response,verbose=verbose,
                                           vel_corr=vel_corr, wlen=wlen,
                                           wfrac=wfrac,
                                           slx=slx, sly=sly, sls=sls)

    def fk_analysis(self, stream, frqlow, frqhigh,
                    prefilter=True, plots=(), static3d=False,correct_3dplane=False,
                    sec_km=False,array_response=False,
                    vel_corr=4.8, wlen=-1, wfrac=0.8,verbose=False,
                    slx=(-10, 10), sly=(-10, 10), sls=0.5):
        """
        FK analysis.

        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :type static3d: bool
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :param array_response: Specify if array response should be used.
        :type array_response: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param wfrac: fraction of sliding window to use for step.
        :type wfrac: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting otions are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        return self._array_analysis_helper(stream=stream, method="FK",
                                           frqlow=frqlow, frqhigh=frqhigh,
                                           prefilter=prefilter, plots=plots,
                                           static3d=static3d,correct_3dplane=correct_3dplane,sec_km=sec_km,
                                           array_r=array_response,vel_corr=vel_corr,
                                           wlen=wlen, wfrac=wfrac,verbose=verbose,
                                           slx=slx, sly=sly, sls=sls)

    def capon_estimator(self, stream, frqlow, frqhigh,
                        prefilter=True, plots=(), static3d=False,correct_3dplane=False,
                        sec_km=False,array_response=False,
                        vel_corr=4.8, wlen=-1, wfrac=0.8,verbose=False,
                        slx=(-10, 10), sly=(-10, 10), sls=0.5):
        """
        Capon's high resolution estimator.

        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :type static3d: bool
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :param array_response: Specify if array response should be used.
        :type array_response: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param wfrac: fraction of sliding window to use for step.
        :type wfrac: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting otions are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        return self._array_analysis_helper(stream=stream, method="CAPON",
                                           frqlow=frqlow, frqhigh=frqhigh,
                                           prefilter=prefilter, plots=plots,
                                           static3d=static3d,correct_3dplane=correct_3dplane,sec_km=sec_km,
                                           array_r=array_response,vel_corr=vel_corr,verbose=verbose,
                                           wlen=wlen, wfrac=wfrac,
                                           slx=slx, sly=sly, sls=sls)

    def _array_analysis_helper(self, stream, method, frqlow, frqhigh,
                               prefilter=True, static3d=False, correct_3dplane=False,sec_km=False,
                               array_r=False, vel_corr=4.8, wlen=-1, wfrac=0.8, slx=(-10, 10),verbose=False,
                               sly=(-10, 10), sls=0.5,
                               plots=()):
        """
        Array analysis wrapper routine.

        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param method: Method used for the array analysis
            (one of "FK": Frequency Wavenumber,
                    "CAPON": Capon's high resolution estimator,
                    "DLS": Delay and Sum,
                    "PWS": Phase Weighted Stack,
                    "SWP": Slowness Whitened Power).
        :type method: str
        :param prefilter: Whether to bandpass data to selected frequency range
        :type prefilter: bool
        :param frqlow: Low corner of frequency range for array analysis
        :type frqlow: float
        :param frqhigh: High corner of frequency range for array analysis
        :type frqhigh: float
        :param static3d: static correction of topography using `vel_corr` as
         velocity (slow!)
        :type static3d: bool
        :param sec_km: switch for input in s/km rather s/deg
        :type sec_km: bool
        :param verbose: Verbose log output
        :type verbose: bool
        :param vel_corr: Correction velocity for static topography correction
         in km/s.
        :type vel_corr: float
        :param wlen: sliding window for analysis in seconds, use -1 to use the
         whole trace without windowing.
        :type wlen: float
        :param wfrac: fraction of sliding window to use for step.
        :type wfrac: float
        :param slx: Min/Max slowness for analysis in x direction [s/km].
        :type slx: (float, float)
        :param sly: Min/Max slowness for analysis in y direction [s/km].
        :type sly: (float, float)
        :param sls: step width of slowness grid [s/km].
        :type sls: float
        :param plots: List or tuple of desired plots that should be plotted for
         each beamforming window.
         Supported options:
         "slowness_baz" for backazimuth-slowness maps for each window,
         "slowness_xy" for slowness_xy maps for each window.
         Further plotting otions are attached to the returned object.
        :rtype: :class:`~obspy.signal.array_analysis.BeamformerResult`
        """
        stream = stream.sort()
        if method not in ("FK", "CAPON", "DLS", "PWS", "SWP"):
            raise ValueError("Invalid method: ''" % method)

        if "slowness_baz" in plots:
            make_slow_map = True
        else:
            make_slow_map = False
        if "slowness_xy" in plots:
            make_slowness_xy = True
        else:
            make_slowness_xy = False

        sllx, slmx = slx
        slly, slmy = sly
        # In terms of a single 'radial' slowness (used e.g. for plotting):
        if sllx < 0 or slly < 0:
            sll = 0
        else:
            sll = np.sqrt(sllx ** 2 + slly ** 2)
        slm = np.sqrt(slmx ** 2 + slmy ** 2)

        # Do not modify the given stream in place.
        st_workon = stream.copy()
        # Trim the stream so all traces are present.
        starttime = max([tr.stats.starttime for tr in st_workon])
        endtime = min([tr.stats.endtime for tr in st_workon])
        st_workon.trim(starttime, endtime)
        st_workon.detrend("linear")
        st_workon.detrend("simple")

        self._attach_coords_to_stream(st_workon)

        if prefilter:
            st_workon.filter('bandpass', freqmin=frqlow, freqmax=frqhigh,
                             zerophase=True)
        else:
            if frqlow is not None or frqhigh is not None:
                warnings.warn("No filtering done. Param 'prefilter' is False.")
        # Making the map plots is efficiently done by saving the power maps to
        # a temporary directory.
        tmpdir = "/tmp/obspy-%s"%UTCDateTime()
        os.mkdir(tmpdir)
        filename_patterns = (os.path.join(tmpdir, 'pow_map_%03d.npy'),
                             os.path.join(tmpdir, 'apow_map_%03d.npy'))
        if make_slow_map or make_slowness_xy:
            def dump(pow_map, apow_map, i):
                np.save(filename_patterns[0] % i, pow_map)
                np.save(filename_patterns[1] % i, apow_map)

        else:
            dump = None

        # Temporarily trim self.inventory so only stations/channels which are
        # actually represented in the traces are kept in the inventory.
        # Otherwise self.geometry and the xyz geometry arrays will have more
        # entries than the stream.
        invbkp = copy.deepcopy(self.inventory)
        self.inventory_cull(st_workon)
        if array_r:
            msll = np.max(np.absolute([sllx, slly, slmx, slmy]))
            frqstep = (frqhigh - frqlow) / 10.
            transff = self.array_transfer_function_freqslowness(msll, sls,
                                                                frqlow,
                                                                frqhigh,
                                                                frqstep)
            print(np.max(transff))
        else:
            transff = None

        try:
            if method == 'FK':
                kwargs = dict(
                    # slowness grid: X min, X max, Y min, Y max, Slow Step
                    sll_x=sllx, slm_x=slmx, sll_y=slly, slm_y=slmy, sl_s=sls,
                    # sliding window properties
                    win_len=wlen, win_frac=wfrac,
                    # frequency properties
                    frqlow=frqlow, frqhigh=frqhigh, prewhiten=0,
                    # restrict output
                    store=dump,
                    semb_thres=-1e9, vel_thres=-1e9, verbose=verbose,
                    # use mlabday to be compatible with matplotlib
                    timestamp='julsec', stime=starttime, etime=endtime,
                    method=0, correct_3dplane=correct_3dplane, vel_cor=vel_corr,
                    static3d=static3d,sec_km=sec_km)

                # here we do the array processing
                start = UTCDateTime()
                outarr = self._covariance_array_processing(st_workon, **kwargs)
                print("Total time in routine: %f\n" % (UTCDateTime() - start))
                t, rel_power, abs_power, baz, slow = outarr.T

            elif method == 'CAPON':
                kwargs = dict(
                    # slowness grid: X min, X max, Y min, Y max, Slow Step
                    sll_x=sllx, slm_x=slmx, sll_y=slly, slm_y=slmy, sl_s=sls,
                    # sliding window properties
                    win_len=wlen, win_frac=wfrac,
                    # frequency properties
                    frqlow=frqlow, frqhigh=frqhigh, prewhiten=0,
                    # restrict output
                    store=dump,
                    semb_thres=-1e9, vel_thres=-1e9, verbose=verbose,
                    # use mlabday to be compatible with matplotlib
                    timestamp='julsec', stime=starttime, etime=endtime,
                    method=1, correct_3dplane=correct_3dplane, vel_cor=vel_corr,
                    static3d=static3d,sec_km=sec_km)

                # here we do the array processing
                start = UTCDateTime()
                outarr = self._covariance_array_processing(st_workon, **kwargs)
                print("Total time in routine: %f\n" % (UTCDateTime() - start))
                t, rel_power, abs_power, baz, slow = outarr.T

            else:
                kwargs = dict(
                    # slowness grid: X min, X max, Y min, Y max, Slow Step
                    sll_x=sllx, slm_x=slmx, sll_y=slly, slm_y=slmy, sl_s=sls,
                    # sliding window properties
                    # frequency properties
                    frqlow=frqlow, frqhigh=frqhigh,
                    # restrict output
                    store=dump,
                    win_len=wlen, win_frac=0.5,
                    nthroot=4, method=method,
                    verbose=verbose, timestamp='julsec',correct_3dplane = correct_3dplane,
                    stime=starttime, etime=endtime, vel_cor=vel_corr,
                    static3d=static3d,sec_km=sec_km)

                # here we do the array processing
                start = UTCDateTime()
                outarr = self._beamforming(st_workon, **kwargs)
                print("Total time in routine: %f\n" % (UTCDateTime() - start))
                t, rel_power, abs_power, baz, slow = outarr.T
                print(outarr.T)
                # abs_power = None

            baz[baz < 0.0] += 360
            if wlen < 0:
                # Need to explicitly specify the timestep of the analysis.
                out = BeamformerResult(inventory=self.inventory,
                                       win_starttimes=t,
                                       slowness_range=np.arange(sll, slm, sls),
                                       max_rel_power=rel_power,
                                       max_abs_power=abs_power,
                                       max_pow_baz=baz, max_pow_slow=slow,
                                       method=method,
                                       timestep=endtime - starttime)
            else:
                out = BeamformerResult(inventory=self.inventory,
                                       win_starttimes=t,
                                       slowness_range=np.arange(sll, slm, sls),
                                       max_rel_power=rel_power,
                                       max_abs_power=abs_power,
                                       max_pow_baz=baz, max_pow_slow=slow,
                                       method=method)

            # now let's do the plotting
            if "slowness_baz" in plots:
                plot_array_analysis(outarr, transff, sllx, slmx, slly, slmy,
                                    sls, filename_patterns, True,
                                    method, array_r, st_workon, starttime,
                                    wlen, endtime,sec_km)
                plt.show()
            if "slowness_xy" in plots:
                plot_array_analysis(outarr, transff, sllx, slmx, slly, slmy,
                                    sls, filename_patterns, False, method,
                                    array_r, st_workon, starttime, wlen,
                                    endtime,sec_km)

                plt.show()
            # Return the beamforming results to allow working more on them,
            # make other plots etc.
            return out
        finally:
            self.inventory = invbkp
            shutil.rmtree(tmpdir)

    def _attach_coords_to_stream(self, stream, origin=None):
        """
        Attaches dictionary with latitude, longitude and elevation to each
        trace in stream as `trace.stats.coords`. Takes into account local
        depth of sensor. If origin is given the distance between the source
        and receiver is calculated and attached as well.

        :param stream: Stream on which the coordiantes are attached to.
        :type stream: :class:`obspy.core.stream.Stream`
        :param origin: Origin of the seismic source.
        :type origin: :class:`~obspy.core.event.origin.Origin`
        """
        geo = self.geometry

        for tr in stream:
            coords = geo[tr.id]
            if origin:
                event_lat = origin.latitude
                event_lng = origin.longitude
                dist = locations2degrees(coords["latitude"],
                                         coords["longitude"], event_lat,
                                         event_lng)
                tr.stats.coordinates = \
                    AttribDict(dict(latitude=coords["latitude"],
                                    longitude=coords["longitude"],
                                    elevation=coords["absolute_height_in_km"],
                                    distance=dist))
            else:
                tr.stats.coordinates = \
                    AttribDict(dict(latitude=coords["latitude"],
                                    longitude=coords["longitude"],
                                    elevation=coords["absolute_height_in_km"]))

    def _covariance_array_processing(self, stream, win_len, win_frac, sll_x,
                                     slm_x, sll_y, slm_y, sl_s, semb_thres,
                                     vel_thres, frqlow, frqhigh, stime, etime,
                                     prewhiten, verbose=False,
                                     timestamp='mlabday', method=0,
                                     correct_3dplane=False, vel_cor=4.,
                                     static3d=False,store=None,sec_km=False):
        """
        Method for FK-Analysis/Capon

        :param stream: Stream object, the trace.stats dict like class must
            contain an :class:`~obspy.core.util.attribdict.AttribDict` with
            'latitude', 'longitude' (in degrees) and 'elevation' (in km), or
            'x', 'y', 'elevation' (in km) items/attributes, as attached in
            `self._array_analysis_helper`/
        :type win_len: float
        :param win_len: Sliding window length in seconds
        :type win_frac: float
        :param win_frac: Fraction of sliding window to use for step
        :type sll_x: float
        :param sll_x: slowness x min (lower)
        :type slm_x: float
        :param slm_x: slowness x max
        :type sll_y: float
        :param sll_y: slowness y min (lower)
        :type slm_y: float
        :param slm_y: slowness y max
        :type sl_s: float
        :param sl_s: slowness step
        :type semb_thres: float
        :param semb_thres: Threshold for semblance
        :type vel_thres: float
        :param vel_thres: Threshold for velocity
        :type frqlow: float
        :param frqlow: lower frequency for fk/capon
        :type frqhigh: float
        :param frqhigh: higher frequency for fk/capon
        :type stime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param stime: Start time of interest
        :type etime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param etime: End time of interest
        :type prewhiten: int
        :param prewhiten: Do prewhitening, values: 1 or 0
        :type timestamp: str
        :param timestamp: valid values: 'julsec' and 'mlabday'; 'julsec'
            returns the timestamp in seconds since 1970-01-01T00:00:00,
            'mlabday' returns the timestamp in days (decimals represent hours,
            minutes and seconds) since '0001-01-01T00:00:00' as needed for
            matplotlib date plotting (see e.g. matplotlib's num2date)
        :type method: int
        :param method: the method to use 0 == bf, 1 == capon
        :param vel_cor: correction velocity (upper layer) in km/s
        :param static3d: a correction of the station height is applied using
            vel_cor the correction is done according to the formula:
            t = rxy*s - rz*cos(inc)/vel_cor
            where inc is defined by inc = asin(vel_cor*slow)
        :type sec_km: bool
        :param sec_km: switch to select betweed s/deg input (false - default) 
               s/km input(true)
        :type store: function
        :param store: A custom function which gets called on each iteration.
            It is called with the relative power map and the time offset as
            first and second arguments and the iteration number as third
            argument. Useful for storing or plotting the map for each
            iteration.
        :return: :class:`numpy.ndarray` of timestamp, relative relpow, absolute
            relpow, backazimuth, slowness
        """
        res = []
        eotr = True

        # check that sampling rates do not vary
        fs = stream[0].stats.sampling_rate
        if len(stream) != len(stream.select(sampling_rate=fs)):
            msg = ('in array-processing sampling rates of traces in stream are'
                   ' not equal')
            raise ValueError(msg)

        grdpts_x = int(((slm_x - sll_x) / sl_s + 0.5) + 1)
        grdpts_y = int(((slm_y - sll_y) / sl_s + 0.5) + 1)

        if correct_3dplane:
            self._correct_with_3dplane(self.geometry)

        if verbose:
            print("geometry:")
            print(self.geometry)
            print("stream contains following traces:")
            print(stream)
            print("stime = " + str(stime) + ", etime = " + str(etime))

        time_shift_table = self._get_timeshift(sll_x, sll_y, sl_s,
                                               grdpts_x, grdpts_y,
                                               vel_cor=vel_cor,
                                               static3d=static3d,sec_km=sec_km)

        spoint, _epoint = _get_stream_offsets(stream, stime, etime)

        # loop with a sliding window over the dat trace array and apply bbfk
        nstat = len(stream)
        fs = stream[0].stats.sampling_rate
        if win_len < 0.:
            nsamp = int((etime - stime) * fs)
            print(nsamp)
            nstep = 1
        else:
            nsamp = int(win_len * fs)
            nstep = int(nsamp * win_frac)

        # generate plan for rfftr
        nfft = next_pow_2(nsamp)
        deltaf = fs / float(nfft)
        nlow = int(frqlow / float(deltaf) + 0.5)
        nhigh = int(frqhigh / float(deltaf) + 0.5)
        nlow = max(1, nlow)  # avoid using the offset
        nhigh = min(nfft // 2 - 1, nhigh)  # avoid using nyquist
        nf = nhigh - nlow + 1  # include upper and lower frequency
        # to speed up the routine a bit we estimate all steering vectors in
        # advance
        steer = np.empty((nf, grdpts_x, grdpts_y, nstat), dtype=np.complex128)
        clibsignal.calcSteer(nstat, grdpts_x, grdpts_y, nf, nlow,
                             deltaf, time_shift_table, steer)
        r = np.empty((nf, nstat, nstat), dtype=np.complex128)
        ft = np.empty((nstat, nf), dtype=np.complex128)
        newstart = stime
        # 0.22 matches 0.2 of historical C bbfk.c
        tap = cosine_taper(nsamp, p=0.22)
        offset = 0
        count = 0  # iteration of loop
        relpow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)
        abspow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)
        while eotr:
            try:
                for i, tr in enumerate(stream):
                    dat = tr.data[spoint[i] + offset:
                                  spoint[i] + offset + nsamp]
                    dat = (dat - dat.mean()) * tap
                    ft[i, :] = np.fft.rfft(dat, nfft)[nlow:nlow + nf]
            except IndexError:
                break
            ft = np.ascontiguousarray(ft, np.complex128)
            relpow_map.fill(0.)
            abspow_map.fill(0.)
            # computing the covariances of the signal at different receivers
            dpow = 0.
            for i in range(nstat):
                for j in range(i, nstat):
                    r[:, i, j] = ft[i, :] * ft[j, :].conj()
                    if method == 1:
                        r[:, i, j] /= np.abs(r[:, i, j].sum())
                    if i != j:
                        r[:, j, i] = r[:, i, j].conjugate()
                    else:
                        dpow += np.abs(r[:, i, j].sum())
            dpow *= nstat
            if method == 1:
                # P(f) = 1/(e.H r(f)^-1 e)
                for n in range(nf):
                    r[n, :, :] = np.linalg.pinv(r[n, :, :], rcond=1e-6)

            errcode = clibsignal.generalizedBeamformer(
                relpow_map, abspow_map, steer, r, nstat, prewhiten,
                grdpts_x, grdpts_y, nf, dpow, method)
            if errcode != 0:
                msg = 'generalizedBeamforming exited with error %d'
                raise Exception(msg % errcode)
            ix, iy = np.unravel_index(relpow_map.argmax(), relpow_map.shape)
            relpow, abspow = relpow_map[ix, iy], abspow_map[ix, iy]
            if store is not None:
                store(relpow_map, abspow_map, count)
            count += 1

            # here we compute baz, slow
            slow_x = sll_x + ix * sl_s
            slow_y = sll_y + iy * sl_s

            slow = np.sqrt(slow_x ** 2 + slow_y ** 2)
            if slow < 1e-8:
                slow = 1e-8
            # Transform from a slowness grid to polar coordinates.
            azimut = 180 * math.atan2(slow_x, slow_y) / math.pi
            baz = azimut % -360 + 180
            if relpow > semb_thres and 1. / slow > vel_thres:
                if timestamp == 'julsec':
                    outtime = newstart
                elif timestamp == 'mlabday':
                    # 719163 == days between 1970 and 0001 + 1
                    outtime = UTCDateTime(newstart.timestamp /
                                          (24. * 3600) + 719163)
                else:
                    msg = "Option timestamp must be one of 'julsec'," \
                          " or 'mlabday'"
                    raise ValueError(msg)
                res.append(np.array([outtime, relpow, abspow, baz,
                                     slow]))
                if verbose:
                    print(newstart, (newstart + (nsamp / fs)), res[-1][1:])
            if (newstart + (nsamp + nstep) / fs) > etime:
                eotr = False
            offset += nstep

            newstart += nstep / fs
        return np.array(res)

    @staticmethod
    def _three_c_dowhiten(fcoeffz, fcoeffn, fcoeffe, deltaf, whiten):
        """
        Amplitude spectra whitening with moving average and window width ww
        and weighting factor: 1/((Z+E+N)/3)
        """
        for nst in range(fcoeffz.shape[0]):
            for nwin in range(fcoeffz.shape[1]):
                ampz = np.abs(fcoeffz[nst, nwin, :])
                ampn = np.abs(fcoeffn[nst, nwin, :])
                ampe = np.abs(fcoeffe[nst, nwin, :])
                # window width can be chosen but must be at least 2 and even:
                ww = int(round(whiten / deltaf))
                if ww == 0:
                    ww = 2
                if ww % 2:
                    ww += 1
                n_freqs = len(ampz)
                csamp = np.zeros((n_freqs, 3), dtype=ampz.dtype)
                csamp[:, 0] = np.cumsum(ampz)
                csamp[:, 1] = np.cumsum(ampe)
                csamp[:, 2] = np.cumsum(ampn)
                ampw = np.zeros(n_freqs, dtype=csamp.dtype)
                for k in range(3):
                    ampw[ww // 2:n_freqs - ww // 2] += \
                        (csamp[ww:, k] - csamp[:-ww, k]) / ww
                # Fill zero elements at start and end of array with closest
                # non-zero value.
                ampw[n_freqs - ww // 2:] = ampw[n_freqs - ww // 2 - 1]
                ampw[:ww // 2] = ampw[ww // 2]
                ampw *= 1 / 3.
                # Weights are 1/ampw unless ampw is very small, then 0.
                weight = np.where(ampw > np.finfo(np.float64).eps * 10.,
                                  1. / (ampw + np.finfo(np.float64).eps), 0.)
                fcoeffz[nst, nwin, :] *= weight
                fcoeffe[nst, nwin, :] *= weight
                fcoeffn[nst, nwin, :] *= weight
        return fcoeffz, fcoeffn, fcoeffe

    def _three_c_do_bf(self, stream_n, stream_e, stream_z, win_len, win_frac,
                       u, sub_freq_range, n_min_stns, polarisation,
                       whiten, phaseonly, coherency, win_average,
                       datalen_sec, uindex, verbose=False):
        """
        Three component beamforming wrapped routine.


        :param n_min_stns: Minimum number of stations for which data must be
         present in a time window, otherwise that window is skipped.
        :param win_frac: fraction of sliding window to use for step


        :return: A :class:`~obspy.signal.array_analysis.BeamformerResult`
        object containing the beamforming results, with dimensions of
        backazimuth range, slowness range, number of windows and number of
        discrete frequencies; as well as frequency and incidence angle arrays
        (the latter will be zero for radial and transversal polarization)

        :param stream_n: Stream of all traces for the North component.
        :param stream_e: Stream of East components.
        :param stream_z: Stream of Up components. Will be ignored for Love
         waves.
        :param win_len: Window length in seconds
        :param win_frac: Overlapping fraction of the window.
        :param u: Array of slowness values.
        :param sub_freq_range: Frequency band (min, max) that is used for
         beamforming and returned. Ideally, use the frequency band of the
         pre-filter.
        :param n_min_stns: Minimum number of stations for which data must be
        :param polarisation: Numeric key specifying the wave polarisation:
                             transvers: 0
                             radial: 1
                             elliptic_retrograde: 2
                             elliptic_prograde: 3
                             p: 4
                             sv' 5
        :param whiten: If set to a number, the 3-component data spectra are
         jointly whitened along the frequency axis with a moving window of
         frequency width 'whiten'.
        :param phaseonly: Whether to totally disregard data amplitudes.
        :param coherency: whether to normalise the beam power spectral density
         by the average station power spectral density of all components
        :param win_average: number of windows to average covariance matrix over
        :param datalen_sec: Difference between start time and end time in [s].
        :param uindex: Slowness range for angle measurments.
        :param verbose: Produce detailed logging information.
        :return: Beamforming results as array, frequencies array,
                 incidence array, window start times array
        """
        # backazimuth range to search
        theo_backazi = np.arange(0, 362, 2) * math.pi / 180.

        # Number of stations should be the same as the number of traces,
        # given the checks in the calling method.
        n_stats = len(stream_n.traces)
        npts = stream_n[0].stats.npts

        geo_array = self._geometry_dict_to_array(
            self._get_geometry_xyz(**self.center_of_gravity))
        # NB at this point these offset arrays will contain three times as many
        # entries as needed because each channel is listed individually. These
        # are cut later on by indexing with the ans array which sorts and
        # selects only the relevant entries.
        x_offsets = geo_array[:, 0]
        y_offsets = geo_array[:, 1]
        # This must be sorted the same as the entries in geo_array!
        # (or channel names, really)
        geo_items_names = []
        for key in sorted(self.geometry):
            geo_items_names.append(key)
        # This is necessary to use np.where below...
        geo_items_names = np.array(geo_items_names)

        # Arrays to hold all traces' data in one:
        _alldata_z = np.zeros((n_stats, npts)) * np.nan
        _alldata_e = _alldata_z.copy()
        _alldata_n = _alldata_z.copy()
        # Array used for sorting and selecting: So far, x_offsets contains
        # offsets for all channels, but the method needs only stations offsets
        # (i.e. a third of the length of the offset array).
        ans = []
        for i, (tr_N, tr_E, tr_Z) in enumerate(zip(stream_n, stream_e,
                                                   stream_z)):
            ans.append(np.where(geo_items_names == tr_N.id)[0][0])
            _alldata_n[i, :] = tr_N.data
            _alldata_e[i, :] = tr_E.data
            _alldata_z[i, :] = tr_Z.data

        fs = stream_n.traces[0].stats.sampling_rate
        # Use np.int to not get the newint type, which causes an error with old
        # versions of numpy (1.6.2 as listed in the minimum requirements).
        nsamp = np.int(win_len * fs)
        # Number of samples to move forward by during a step.
        nstep = int(nsamp * win_frac)
        # Number of windows is determined by data length minus one window
        # length divided by step length, then adding the one omitted window.
        num_win = int((datalen_sec * fs - nsamp) / nstep) + 1
        out_wins = int(np.floor(num_win / win_average))
        if not out_wins > 0:
            msg = "Zero output windows! Check data length, and parameters " \
                  "window_length, win_frac and win_average."
            raise ValueError(msg)

        alldata_z = np.zeros((n_stats, num_win, nsamp))
        alldata_n, alldata_e = alldata_z.copy(), alldata_z.copy()
        nst = np.zeros(num_win)

        # Iterate over the beamfoming windows:
        for i in range(num_win):
            for n in range(n_stats):
                if not np.isnan(_alldata_z[n, i * nstep:i * nstep +
                                           nsamp]).any() \
                        and not np.isnan(_alldata_n[n, i * nstep:i * nstep +
                                                    nsamp]).any() \
                        and not np.isnan(_alldata_e[n, i * nstep:i * nstep +
                                                    nsamp]).any():
                    # All data, tapered.
                    alldata_z[n, i, :] = _alldata_z[n, i * nstep:
                                                    i * nstep + nsamp] * \
                        cosine_taper(nsamp)
                    alldata_n[n, i, :] = _alldata_n[n, i * nstep:
                                                    i * nstep + nsamp] * \
                        cosine_taper(nsamp)
                    alldata_e[n, i, :] = _alldata_e[n, i * nstep:
                                                    i * nstep + nsamp] * \
                        cosine_taper(nsamp)
                    nst[i] += 1

        # Need an array of the starting times of the beamforming windows for
        # later reference (e.g. plotting). The precision of these should not
        # exceed what is reasonable given the sampling rate: the different
        # traces have the same start times only to within a sample.
        avg_starttime = UTCDateTime(np.mean([tr.stats.starttime.timestamp for
                                             tr in stream_n.traces]))
        window_start_times = \
            np.array([UTCDateTime(avg_starttime + i * nstep / fs,
                                  precision=len(str(fs).split('.')[0]))
                      for i in range(num_win) if i % win_average == 0])
        if verbose:
            print(nst, ' stations/window; average over ', win_average)

        # Do Fourier transform.
        deltat = stream_n.traces[0].stats.delta
        freq_range = np.fft.fftfreq(nsamp, deltat)
        # Use a narrower 'frequency range' of interest for evaluating incidence
        # angle.
        lowcorner = sub_freq_range[0]
        highcorner = sub_freq_range[1]
        index = np.where((freq_range >= lowcorner) &
                         (freq_range <= highcorner))[0]
        fr = freq_range[index]

        # for final Power Spectral Density output using half-sided spectrum,
        # traces are normalized by SQRT(fs*windowing-function factor*0.5)
        fcoeffz = np.fft.fft(alldata_z, n=nsamp, axis=-1) / \
            np.sqrt(fs * (cosine_taper(nsamp) ** 2).sum() * 0.5)
        fcoeffn = np.fft.fft(alldata_n, n=nsamp, axis=-1) / \
            np.sqrt(fs * (cosine_taper(nsamp) ** 2).sum() * 0.5)
        fcoeffe = np.fft.fft(alldata_e, n=nsamp, axis=-1) / \
            np.sqrt(fs * (cosine_taper(nsamp) ** 2).sum() * 0.5)
        fcoeffz = fcoeffz[:, :, index]
        fcoeffn = fcoeffn[:, :, index]
        fcoeffe = fcoeffe[:, :, index]
        deltaf = 1. / (nsamp * deltat)

        if whiten:
            try:
                float(whiten)
            except ValueError or TypeError:
                msg = 'Whiten parameter must be digits. It was set to 0.01.'
                warnings.warn(msg)
                whiten = 0.01
            if whiten >= fr[-1] - fr[0]:
                msg = ('Moving frequency window is %s, it equals or exceeds '
                       'the entire frequency range and was set to 0.01 now.')
                warnings.warn(msg % whiten)
                whiten = 0.01
            fcoeffz, fcoeffn, fcoeffe = self._three_c_dowhiten(
                fcoeffz, fcoeffn, fcoeffe, deltaf, whiten)
        if phaseonly:
            fcoeffz = np.exp(1j * np.angle(fcoeffz))
            fcoeffn = np.exp(1j * np.angle(fcoeffn))
            fcoeffe = np.exp(1j * np.angle(fcoeffe))

        # slowness vector u and slowness vector component scale u_x and u_y
        theo_backazi = theo_backazi.reshape((theo_backazi.size, 1))
        u_y = -np.cos(theo_backazi)
        u_x = -np.sin(theo_backazi)

        # vector of source direction dependent plane wave travel-distance to
        # reference point (positive value for later arrival/negative for
        # earlier arr)
        x_offsets = np.array(x_offsets)
        y_offsets = np.array(y_offsets)
        # This sorts the offset value arrays.
        x_offsets = x_offsets[np.array(ans)]
        y_offsets = y_offsets[np.array(ans)]
        # The steering vector corresponds to the 'rho_m * cos(theta_m - theta)'
        # factor in Table 1 of Esmersoy et al., 1985.
        steering = u_y * y_offsets + u_x * x_offsets

        # polarizations [Z,E,N]
        # incident angle or atan(H/V)
        incs = np.arange(5, 90, 10) * math.pi / 180.

        def pol_transverse(azi):
            pol_e = math.cos(theo_backazi[azi])
            pol_n = -1. * math.sin(theo_backazi[azi])
            return pol_e, pol_n

        def pol_rayleigh_retro(azi):
            pol_e = math.sin(theo_backazi[azi])
            pol_n = math.cos(theo_backazi[azi])
            return pol_e, pol_n

        def pol_rayleigh_prog(azi):
            pol_e = -1 * math.sin(theo_backazi[azi])
            pol_n = -1 * math.cos(theo_backazi[azi])
            return pol_e, pol_n

        def pol_p(azi):
            pol_e = -1 * math.sin(theo_backazi[azi])
            pol_n = -1 * math.cos(theo_backazi[azi])
            return pol_e, pol_n

        def pol_sv(azi):
            pol_e = math.sin(theo_backazi[azi])
            pol_n = math.cos(theo_backazi[azi])
            return pol_e, pol_n

        cz = [0., 0., 1j, 1j, 1., 1.]
        ch = (pol_transverse, pol_rayleigh_retro, pol_rayleigh_retro,
              pol_rayleigh_prog, pol_p, pol_sv)

        nfreq = len(fr)
        beamres = np.zeros((len(theo_backazi), u.size,
                            max(out_wins, len(window_start_times)), nfreq))
        incidence = np.zeros((max(out_wins, len(window_start_times)), nfreq))
        win_average = int(win_average)
        for f in range(nfreq):
            omega = 2 * math.pi * fr[f]
            for win in range(0, out_wins * win_average, win_average):
                if any(nst[win:win + win_average] < n_min_stns) or any(
                                nst[win:win + win_average] != nst[win]):
                    continue
                sz = np.squeeze(fcoeffz[:, win, f])
                sn = np.squeeze(fcoeffn[:, win, f])
                se = np.squeeze(fcoeffe[:, win, f])

                y = np.concatenate((sz, sn, se))
                y = y.reshape(1, y.size)
                yt = y.T.copy()
                r = np.dot(yt, np.conjugate(y))

                for wi in range(1, win_average):
                    sz = np.squeeze(fcoeffz[:, win + wi, f])
                    sn = np.squeeze(fcoeffn[:, win + wi, f])
                    se = np.squeeze(fcoeffe[:, win + wi, f])

                    y = np.concatenate((sz, sn, se))
                    y = y.reshape(1, y.size)
                    yt = y.T.copy()
                    r += np.dot(yt, np.conjugate(y))

                r /= float(win_average)

                res = np.zeros((len(theo_backazi), len(u), len(incs)))
                for vel in range(len(u)):
                    e_steer = np.exp(-1j * steering * omega * u[vel])
                    e_steere = e_steer.copy()
                    e_steern = e_steer.copy()
                    e_steere = (e_steere.T * np.array([ch[polarisation](azi)[0]
                                for azi in range(len(theo_backazi))])).T
                    e_steern = (e_steern.T * np.array([ch[polarisation](azi)[1]
                                for azi in range(len(theo_backazi))])).T

                    if polarisation in [0, 1]:
                        w = np.concatenate(
                            (e_steer * cz[polarisation], e_steern, e_steere),
                            axis=1)
                        wt = w.T.copy()
                        beamres[:, vel, int(win / win_average), f] = 1. / (
                                nst[win] * nst[win]) * abs(
                                (np.conjugate(w) * np.dot(r, wt).T).sum(1))
                        if coherency:
                            beamres[:, vel, int(win / win_average), f] /= \
                                abs(np.sum(np.diag(r)))

                    elif polarisation in [2, 3, 4]:
                        for inc_angle in range(len(incs)):
                            w = np.concatenate((e_steer * cz[polarisation] *
                                                np.cos(incs[inc_angle]),
                                                e_steern *
                                                np.sin(incs[inc_angle]),
                                                e_steere *
                                                np.sin(incs[inc_angle])),
                                               axis=1)
                            wt = w.T.copy()
                            res[:, vel, inc_angle] = 1. / (
                                    nst[win] * nst[win]) * abs(
                                    (np.conjugate(w) * np.dot(r, wt).T).sum(1))
                            if coherency:
                                res[:, vel, inc_angle] /= \
                                    abs(np.sum(np.diag(r)))

                    elif polarisation == 5:
                        for inc_angle in range(len(incs)):
                            w = np.concatenate((e_steer * cz[polarisation] *
                                                np.sin(incs[inc_angle]),
                                                e_steern *
                                                np.cos(incs[inc_angle]),
                                                e_steere *
                                                np.cos(incs[inc_angle])),
                                               axis=1)
                            wt = w.T.copy()
                            res[:, vel, inc_angle] = 1. / (
                                    nst[win] * nst[win]) * abs(
                                    (np.conjugate(w) * np.dot(r, wt).T).sum(1))
                            if coherency:
                                res[:, vel, inc_angle] /= \
                                    abs(np.sum(np.diag(r)))

                if polarisation > 1:
                    i, j, k = np.unravel_index(np.argmax(res[:, uindex, :]),
                                               res.shape)
                    beamres[:, :, int(win / win_average), f] = res[:, :, k]
                    incidence[int(win / win_average), f] = incs[k] * 180. / \
                        math.pi

        return beamres, fr, incidence, window_start_times

    def three_component_beamforming(self, stream_n, stream_e, stream_z, wlen,
                                    smin, smax, sstep, wavetype, freq_range,
                                    n_min_stns=5, win_average=1, win_frac=1,
                                    whiten=False, phaseonly=False,
                                    coherency=False):
        """
        Do three-component beamforming following [Esmersoy1985]_.

        Three streams representing N, E, Z oriented components must be given,
        where the traces contained are from the different stations. The
        traces must all have same length and start/end times (to within
        sampling distance). (hint: check length with trace.stats.npts)
        The given streams are not modified in place. All trimming, filtering,
        downsampling should be done previously.
        The beamforming can distinguish horizontally transversal (SH), radial,
        prograde/retrograde elliptical, longitudinal (P) and vertically
        transversal (SV) polarization and performs grid searchs over slowness,
        azimuth and incidence angle, respectively arctangent of the H/V ratio.
        Station location information is taken from the array's inventory, so
        that must contain station or channel location information about all
        traces used (or more, the inventory is then non-permanently 'pruned').
        NB all channels of a station must be located in the same location for
        this method.

        :param stream_n: Stream of all traces for the North component.
        :param stream_e: Stream of East components.
        :param stream_z: Stream of Up components. Will be ignored for Love
         waves.
        :param wlen: window length in seconds
        :param smin: minimum slowness of the slowness grid [s/km]
        :param smax: maximum slowness [s/km]
        :param sstep: slowness step [s/km]
        :param wavetype: 'transvers', 'radial', 'elliptic_retrograde',
         'elliptic_prograde', 'P', or 'SV'
        :param freq_range: Frequency band (min, max) that is used for
         beamforming and returned. Ideally, use the frequency band of the
         pre-filter.
        :param n_min_stns: Minimum number of stations for which data must be
         present in a time window, otherwise that window is skipped.
        :param win_average: number of windows to average covariance matrix over
        :param win_frac: fraction of sliding window to use for step
        :param whiten: if set to a number, the 3-component data spectra are
         jointly whitened along the frequency axis with a moving window of
         frequency width 'whiten'
        :param phaseonly: whether to totally disregard data amplitudes
        :param coherency: whether to normalise the beam power spectral density
         by the average station power spectral density of all components
        :return: A :class:`~obspy.signal.array_analysis.BeamformerResult`
        object containing the beamforming results, with dimensions of
        backazimuth range, slowness range, number of windows and number of
        discrete frequencies; as well as frequency and incidence angle arrays
        (the latter will be zero for radial and transversal polarization).
        """
        pol_dict = {'transvers': 0, 'radial': 1, 'elliptic_retrograde': 2,
                    'elliptic_prograde': 3, 'p': 4, 'sv': 5}
        if wavetype.lower() not in pol_dict:
            raise ValueError('Invalid option for wavetype: {}'
                             .format(wavetype))
        if len(set(len(vel.traces) for vel in (stream_n, stream_e,
                                               stream_z))) > 1:
            raise ValueError("All three streams must have same number of "
                             "traces.")
        if len(stream_n.traces) == 0:
            raise ValueError("Streams do not seem to contain any traces.")

        # from _array_analysis_helper:
        starttime = max(max([tr.stats.starttime for tr in st]) for st in
                        (stream_n, stream_e, stream_e))
        min_starttime = min(min([tr.stats.starttime for tr in st]) for st in
                            (stream_n, stream_e, stream_e))
        endtime = min(min([tr.stats.endtime for tr in st]) for st in
                      (stream_n, stream_e, stream_e))
        max_endtime = max(max([tr.stats.endtime for tr in st]) for st in
                          (stream_n, stream_e, stream_e))

        delta_common = stream_n.traces[0].stats.delta
        npts_common = stream_n.traces[0].stats.npts
        if max(abs(min_starttime - starttime),
               abs(max_endtime - endtime)) > delta_common:
            raise ValueError("Traces do not have identical start/end times. "
                             "Trim to same times (within sample accuracy) "
                             "and ensure all traces have the same number "
                             "of samples!")

        # Check for equal deltas and number of samples:
        for st in (stream_n, stream_e, stream_z):
            for tr in st:
                if tr.stats.npts != npts_common:
                    raise ValueError('Traces do not have identical number of '
                                     'samples.')
                if tr.stats.delta != delta_common:
                    raise ValueError('Traces do not have identical sampling '
                                     'rates.')
        datalen_sec = endtime - starttime

        # Sort all traces just to make sure they're in the same order.
        for st in (stream_n, stream_e, stream_z):
            st.sort()

        for trN, trE, trZ in zip(stream_n, stream_e, stream_z):
            if len(set('{}.{}'.format(tr.stats.network, tr.stats.station)
                       for tr in (trN, trE, trZ))) > 1:
                raise ValueError("Traces are not from same stations.")

        # Temporarily trim self.inventory so only stations/channels which are
        # actually represented in the traces are kept in the inventory.
        # Otherwise self.geometry and the xyz geometry arrays will have more
        # entries than the stream.
        invbkp = copy.deepcopy(self.inventory)
        allstreams = stream_n + stream_e + stream_z
        self.inventory_cull(allstreams)

        try:
            if wlen < smax * self.aperture:
                raise ValueError('Window length is smaller than maximum given'
                                 ' slowness times aperture.')
            # s/km  slowness range calculated
            u = np.arange(smin, smax, sstep)
            # Slowness range evaluated for (incidence) angle measurement
            # (Rayleigh, P, SV):
            # These values are a bit arbitrary for now:
            uindex = np.where((u > 0.5 * smax + smin) &
                              (u < 0.8 * smax + smin))[0]

            bf_results, freqs, incidence, window_start_times = \
                self._three_c_do_bf(stream_n, stream_e, stream_z,
                                    win_len=wlen, win_frac=win_frac, u=u,
                                    sub_freq_range=freq_range,
                                    n_min_stns=n_min_stns,
                                    polarisation=pol_dict[wavetype.lower()],
                                    whiten=whiten,
                                    phaseonly=phaseonly,
                                    coherency=coherency,
                                    win_average=win_average,
                                    datalen_sec=datalen_sec,
                                    uindex=uindex)

            out = BeamformerResult(inventory=self.inventory,
                                   win_starttimes=window_start_times,
                                   slowness_range=u, full_beamres=bf_results,
                                   freqs=freqs, incidence=incidence,
                                   method='3C ({})'.format(wavetype),
                                   timestep=wlen * win_frac * win_average)

        finally:
            self.inventory = invbkp

        return out

    def plot_radial_transfer_function(self, smin, smax, sstep, freqs):
        """
        Plot array transfer function radially, as function of slowness.

        :param smin: Minimum slowness value.
        :param smax: Maximum slowness value.
        :param sstep: Slowness step.
        :param freqs: List of frequencies for which the radial transfer
                      function should be plotted.
        """
        u = np.arange(smin, smax, sstep)
        theo_backazi = np.arange(0, 362, 2) * math.pi / 180.
        theo_backazi = theo_backazi.reshape((theo_backazi.size, 1))
        u_y = -np.cos(theo_backazi)
        u_x = -np.sin(theo_backazi)
        geo_array = self._geometry_dict_to_array(
            self._get_geometry_xyz(**self.center_of_gravity))
        x_ = geo_array[:, 0]
        y_ = geo_array[:, 1]
        x_ = np.array(x_)
        y_ = np.array(y_)
        steering = u_y * y_ + u_x * x_
        theo_backazi = theo_backazi[:, 0]
        beamres = np.zeros((len(theo_backazi), u.size))
        for f in freqs:
            omega = 2. * math.pi * f
            r = np.ones((steering.shape[1], steering.shape[1]))
            for vel in range(len(u)):
                w = np.exp(-1j * steering * omega * u[vel])
                wt = w.T.copy()
                beamres[:, vel] = 1. / (
                    steering.shape[1] * steering.shape[1]) * abs(
                    (np.conjugate(w) * np.dot(r, wt).T).sum(1))
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection='polar')
            cmap = plt.cm.get_cmap('viridis')
            contf = ax.contourf(theo_backazi, u,
                                beamres.T, 40, cmap=cmap, antialiased=True)
            ax.contour(theo_backazi, u,
                       beamres.T, 40, cmap=cmap)
            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            ax.set_rmax(u[-1])
            # This means that if u does not start at 0, the plot will show a
            # hole in the middle rather than stitching it up.
            ax.set_rmin(-0)
            fig.colorbar(contf)
            ax.grid(True)
            ax.set_title('Transfer function at f = ' + str(f) + '.')
            plt.tight_layout()
        plt.show()

    def plot_transfer_function_wavenumber(self, klim, kstep):
        """
        Plot array transfer function as function of wavenumber.

        :param klim: Maximum wavenumber (symmetric about zero in
         x and y directions).
        :param kstep: Step in wavenumber.
        """
        transff = self.array_transfer_function_wavenumber(klim, kstep)
        self._plot_transfer_function_helper(transff, klim, kstep)
        plt.xlabel('Wavenumber West-East')
        plt.ylabel('Wavenumber North-South')
        plt.show()

    def plot_transfer_function_freqslowness(self, slim, sstep, freq_min,
                                            freq_max, freq_step):
        """
        Plot array transfer function as function of slowness and frequency.

        :param slim: Maximum slowness (symmetric about zero in x and y
         directions).
        :param sstep: Step in slowness.
        :param freq_min: Minimum frequency in signal.
        :param freq_max: Maximum frequency in signal.
        :param freq_step: Frequency sample distance.
        """
        transff = self.array_transfer_function_freqslowness(
            slim, sstep, freq_min, freq_max, freq_step)
        self._plot_transfer_function_helper(transff, slim, sstep)
        plt.xlabel('Slowness West-East')
        plt.ylabel('Slowness North-South')
        plt.show()

    @staticmethod
    def _plot_transfer_function_helper(transff, lim, step):
        """
        Plot array transfer function.

        :param transff: Transfer function to plot.
        :param lim: Maximum value of slowness/wavenumber.
        :param step: Step in slowness/wavenumber.
        """
        ranges = np.arange(-lim, lim + step, step)
        # plt.pcolor(ranges, ranges, transff.T, cmap=cm.viridis)
        plt.contour(ranges, ranges, transff.T, 10)
        plt.colorbar()
        plt.clim(vmin=0., vmax=1.)
        plt.xlim(-lim, lim)
        plt.ylim(-lim, lim)

    def array_transfer_function_wavenumber(self, klim, kstep):
        """
        Return array transfer function as a function of wavenumber difference.

        :param klim: Either a float to use symmetric limits for wavenumber
            differences or the tuple (kxmin, kxmax, kymin, kymax).
        :param kstep: Step in wavenumber.
        """
        return self._array_transfer_function_helper(klim, kstep, 'wavenumber')

    def array_transfer_function_freqslowness(self, slim, sstep, fmin, fmax,
                                             fstep):
        """
        Return array transfer function as a function of slowness difference
        and frequency.

        :param slim: Either a float to use symmetric limits for slowness
            differences or the tuple (sxmin, sxmax, symin, symax).
        :param sstep: Step in frequency.
        :param fmin: Minimum frequency in signal.
        :param fmax: Maximum frequency in signal.
        :param fstep: Frequency sample distance.
        """

        return self._array_transfer_function_helper(
            slim, sstep, 'slowness', fmin, fmax, fstep)

    def _array_transfer_function_helper(
            self, plim, pstep, param, fmin=None, fmax=None, fstep=None):
        """
        Return array transfer function as function of wavenumber or slowness
        and frequency.

        :param plim: Either a float to use symmetric limits for slowness/
            wavenumber differences or the tuple (pxmin, pxmax, sxmin, sxmax).
        :param pstep: Step in wavenumber/slowness.
        :param param: 'wavenumber' or 'slowness'
        :param fmin: Minimum frequency (only for slowness calculation).
        :param fmax: Maximum frequency (only for slowness calculation).
        :param fstep: Frequency sample distance (only with slowness).
        :return: Array transfer function as function of either wavenumber or
                 slowness and frequency
        """
        if isinstance(plim, (float, int, np.int32, np.int64)):
            pxmin = -plim
            pxmax = plim
            pymin = -plim
            pymax = plim
        elif len(plim) == 4:
            pxmin = plim[0]
            pxmax = plim[1]
            pymin = plim[2]
            pymax = plim[3]
        else:
            raise TypeError('Parameter slim must either be a float '
                            'or a tuple of length 4.')
        # geometry = self._geometry_dict_to_array(self._get_geometry_xyz(
        #    **self.center_of_gravity))
        geometry = self._geometry_dict_to_array(self.geometry)
        npx = int(np.ceil((pxmax + pstep / 10. - pxmin) / pstep))
        npy = int(np.ceil((pymax + pstep / 10. - pymin) / pstep))
        transff = np.empty((npx, npy))

        if param == 'wavenumber':
            for i, kx in enumerate(np.arange(pxmin, pxmax + pstep / 10.,
                                             pstep)):
                for j, ky in enumerate(np.arange(pymin, pymax + pstep / 10.,
                                                 pstep)):
                    _sum = 0j
                    for k in range(len(geometry)):
                        _sum += np.exp(complex(0., geometry[k, 0] * kx +
                                               geometry[k, 1] * ky))
                    transff[i, j] = abs(_sum) ** 2

        elif param == 'slowness':
            nf = int(np.ceil((fmax + fstep / 10. - fmin) / fstep))
            buff = np.zeros(nf)
            for i, sx in enumerate(np.arange(pxmin, pxmax + pstep / 10.,
                                             pstep)):
                for j, sy in enumerate(np.arange(pymin, pymax + pstep / 10.,
                                                 pstep)):
                    for m, f in enumerate(np.arange(fmin, fmax + fstep / 10.,
                                                    fstep)):
                        _sum = 0j
                        for n in np.arange(len(geometry)):
                            _sum += np.exp(complex(0., (geometry[n, 0] * sx +
                                                        geometry[n, 1] * sy) *
                                                   2 * np.pi * f))
                        buff[m] = abs(_sum) ** 2
                    transff[i, j] = cumtrapz(buff, dx=fstep)[-1]

        transff /= transff.max()
        return transff

    def _beamforming(self, stream, sll_x, slm_x, sll_y, slm_y, sl_s, frqlow,
                     frqhigh, stime, etime, win_len=-1, win_frac=0.5,
                     verbose=False, timestamp='mlabday',
                     method="DLS", nthroot=1, store=None,
                     correct_3dplane=False, static3d=False, sec_km=False,vel_cor=4.):
        """
        Method for Delay and Sum/Phase Weighted Stack/Whitened Slowness Power

        :param stream: Stream object.
        :param sll_x: slowness x min (lower)
        :param slm_x: slowness x max
        :param sll_y: slowness y min (lower)
        :param slm_y: slowness y max
        :param sl_s: slowness step
        :type stime: UTCDateTime
        :param stime: Starttime of interest
        :type etime: UTCDateTime
        :param etime: Endtime of interest
        :param win_len: length for sliding window analysis, default is -1
         which means the whole trace;
        :param win_frac of win_len which is used to 'hop' forward in time
        :param timestamp: valid values: 'julsec' and 'mlabday'; 'julsec'
         returns the timestamp in secons since 1970-01-01T00:00:00,
         'mlabday' returns the timestamp in days (decimals represent hours,
         minutes and seconds) since '0001-01-01T00:00:00' as needed for
         matplotlib date plotting (see e.g. matplotlibs num2date).
        :param method: the method to use "DLS" delay and sum; "PWS" phase
         weighted stack; "SWP" slowness weightend power spectrum
        :param nthroot: nth-root processing; nth gives the root (1,2,3,4),
         default 1 (no nth-root)
        :type store: function
        :param store: A custom function which gets called on each iteration. It
         is called with the relative power map and the time offset as first and
         second arguments and the iteration number as third argument. Useful
         for storing or plotting the map for each iteration.
        :param correct_3dplane: if Yes than a best (LSQ) plane will be fitted
         into the array geometry. Mainly used with small apature arrays at
         steep flanks.
        :param static3d: if yes the station height of am array station is
         taken into account according to the formula:
            tj = -xj*sxj - yj*syj + zj*cos(inc)/vel_cor
         the inc angle is slowness dependend and thus must
         be estimated for each grid-point:
            inc = asin(v_cor*slow)
        :type sec_km: bool
        :param sec_km: switch for either s/deg (default - false) or s/km input (true
        :param vel_cor: Velocity for the upper layer (static correction)
         in km/s.
        :return: numpy.ndarray of timestamp, relative relpow, absolute relpow,
         backazimut, slowness, maximum beam (for DLS)
        """
        res = []
        eotr = True

        # check that sampling rates do not vary
        fs = stream[0].stats.sampling_rate
        nstat = len(stream)
        if len(stream) != len(stream.select(sampling_rate=fs)):
            msg = 'sampling rates of traces in stream are not equal'
            raise ValueError(msg)

        # loop with a sliding window over the dat trace array and apply bbfk

        grdpts_x = int(((slm_x - sll_x) / sl_s + 0.5) + 1)
        grdpts_y = int(((slm_y - sll_y) / sl_s + 0.5) + 1)

        abspow_map = np.empty((grdpts_x, grdpts_y), dtype='f8')
        geometry = self._geometry_dict_to_array(self._get_geometry_xyz(
            correct_3dplane=correct_3dplane,
            **self.center_of_gravity))

        if verbose:
            print("geometry:")
            print(geometry)
            print("stream contains following traces:")
            print(stream)
            print("stime = " + str(stime) + ", etime = " + str(etime))

        time_shift_table = self._get_timeshift(sll_x, sll_y, sl_s,
                                               grdpts_x, grdpts_y,
                                               vel_cor=vel_cor,
                                               static3d=static3d,sec_km=sec_km)

        mini = np.min(time_shift_table[:, :, :])
        maxi = np.max(time_shift_table[:, :, :])
        spoint, _epoint = _get_stream_offsets(stream, (stime - mini),
                                              (etime - maxi))

        # recalculate the maximum possible trace length
        #    ndat = int(((etime-maxi) - (stime-mini))*fs)
        if win_len < 0:
            nsamp = int(((etime - maxi) - (stime - mini)) * fs)
        else:
            # nsamp = int((win_len-np.abs(maxi)-np.abs(mini)) * fs)
            nsamp = int(win_len * fs)

        if nsamp <= 0:
            print('Data window too small for slowness grid')
            print('Must exit')
            sys.exit()

        nstep = int(nsamp * win_frac)

        stream.detrend()
        newstart = stime
        offset = 0
        count = 0
        while eotr:
            max_beam = 0.
            if method == 'DLS':
                for x in range(grdpts_x):
                    for y in range(grdpts_y):
                        singlet = 0.
                        beam = np.zeros(nsamp, dtype='f8')
                        for i in range(nstat):
                            s = spoint[i] + int(
                                time_shift_table[i, x, y] * fs + 0.5)
                            try:
                                shifted = stream[i].data[s + offset:
                                                         s + nsamp + offset]
                                if len(shifted) < nsamp:
                                    shifted = np.pad(
                                        shifted, (0, nsamp - len(shifted)),
                                        'mean') #'constant, constant_values=(0, 1))
                                singlet += 1. / nstat * np.sum(shifted *
                                                               shifted)
                                beam += 1. / nstat * np.power(
                                    np.abs(shifted), 1. / nthroot) * \
                                    shifted / np.abs(shifted)
                            except IndexError:
                                break
                        beam = np.power(np.abs(beam), nthroot) * \
                            beam / np.abs(beam)
                        beam=np.nan_to_num(beam)
                        bs = np.sum(beam * beam)
                        abspow_map[x, y] = bs / singlet
                        if abspow_map[x, y] > max_beam:
                            max_beam = abspow_map[x, y]
                            beam_max = beam
            if method == 'PWS':
                for x in range(grdpts_x):
                    for y in range(grdpts_y):
                        singlet = 0.
                        beam = np.zeros(nsamp, dtype='f8')
                        stack = np.zeros(nsamp, dtype='c8')
                        for i in range(nstat):
                            s = spoint[i] + int(time_shift_table[i, x, y] *
                                                fs + 0.5)
                            try:
                                shifted = sp.signal.hilbert(stream[i].data[
                                    s + offset: s + nsamp + offset])
                                if len(shifted) < nsamp:
                                    shifted = np.pad(
                                        shifted, (0, nsamp - len(shifted)),
                                        'constant', constant_values=(0, 1))
                            except IndexError:
                                break
                            phase = np.arctan2(shifted.imag, shifted.real)
                            stack.real += np.cos(phase)
                            stack.imag += np.sin(phase)
                        coh = 1. / nstat * np.abs(stack)
                        for i in range(nstat):
                            s = spoint[i] + int(
                                time_shift_table[i, x, y] * fs + 0.5)
                            shifted = stream[i].data[
                                s + offset: s + nsamp + offset]
                            singlet += 1. / nstat * np.sum(shifted * shifted)
                            try:
                                beam += 1. / nstat * shifted * \
                                        np.power(coh, nthroot)
                            except ValueError:
                                _coh = coh[:len(shifted)]
                                temp = 1. / nstat * shifted * \
                                    np.power(_coh, nthroot)
                                beam += np.pad(temp, (0, beam.shape[0] -
                                                      temp.shape[0]),
                                               'constant')
                        beam=np.nan_to_num(beam)
                        bs = np.sum(beam * beam)
                        abspow_map[x, y] = bs / singlet
                        if abspow_map[x, y] > max_beam:
                            max_beam = abspow_map[x, y]
                            beam_max = beam
            if method == 'SWP':
                # generate plan for rfftr
                nfft = next_pow_2(nsamp)
                deltaf = fs / float(nfft)
                nlow = int(frqlow / float(deltaf) + 0.5)
                nhigh = int(frqhigh / float(deltaf) + 0.5)
                nlow = max(1, nlow)  # avoid using the offset
                nhigh = min(nfft / 2 - 1, nhigh)  # avoid using nyquist
                nf = nhigh - nlow + 1  # include upper and lower frequency

                steer = np.empty((int(nf), grdpts_x, grdpts_y, nstat),
                                 dtype='c16')
                spec = np.zeros((nstat, int(nf)), dtype='c16')
                time_shift_table *= -1.
                clibsignal.calcSteer(nstat, grdpts_x, grdpts_y, int(nf), nlow,
                                     deltaf, time_shift_table, steer)
                try:
                    for i in range(nstat):
                        dat = stream[i].data[spoint[i] + offset:
                                             spoint[i] + offset + nsamp]

                        # TODO: cosine_taper(len(dat), p=0.22) or
                        #   cosine_taper(nsamp, p=0.22)[0:len(dat)]
                        tap = cosine_taper(len(dat), p=0.22)
                        dat = (dat - dat.mean()) * tap
                        spec[i, :] = np.fft.rfft(dat, nfft)[int(nlow):
                                                            int(nlow + nf)]
                except IndexError:
                    break

                for i in range(grdpts_x):
                    for j in range(grdpts_y):
                        for m in range(int(nf)):
                            for n in range(nstat):
                                steer[m, i, j, n] *= spec[n, m]

                beam = np.absolute(np.sum(steer, axis=3))
                less = np.max(beam, axis=1)
                max_buffer = np.max(less, axis=1)

                for i in range(grdpts_x):
                    for j in range(grdpts_y):
                        abspow_map[i, j] = np.sum(beam[:, i, j] /
                                                  max_buffer[:],
                                                  axis=0) / float(nf)

                beam_max = stream[0].data[spoint[0] + offset:
                                          spoint[0] + nsamp + offset]

            ix, iy = np.unravel_index(abspow_map.argmax(), abspow_map.shape)
            abspow = abspow_map[ix, iy]
            if store is not None:
                store(abspow_map, beam_max, count)
            count += 1
            # here we compute baz, slow
            slow_x = sll_x + ix * sl_s
            slow_y = sll_y + iy * sl_s

            slow = np.sqrt(slow_x ** 2 + slow_y ** 2)
            if slow < 1e-8:
                slow = 1e-8
            azimut = 180 * math.atan2(slow_x, slow_y) / math.pi
            baz = azimut % -360 + 180
            if timestamp == 'julsec':
                outtime = newstart
            elif timestamp == 'mlabday':
                # 719163 == days between 1970 and 0001 + 1
                outtime = UTCDateTime(newstart.timestamp /
                                      (24. * 3600) + 719163)
            else:
                msg = "Option timestamp must be one of 'julsec'," \
                      " or 'mlabday'"
                raise ValueError(msg)
            res.append(np.array([outtime, abspow, abspow, baz,
                                     slow]))

            if verbose:
                print(newstart, (newstart + (nsamp / fs)), res[-1][1:])
            if (newstart + (nsamp + nstep) / fs) > etime:
                eotr = False
            offset += nstep

            newstart += nstep / fs

        #if timestamp == 'julsec':
        #    pass
        #elif timestamp == 'mlabday':
        #    # 719162 == hours between 1970 and 0001
        #    res[:, 0] = res[:, 0] / (24. * 3600) + 719162
        #else:
        #    msg = "Option timestamp must be one of 'julsec', or 'mlabday'"
        #    raise ValueError(msg)
        return np.array(res)

    @staticmethod
    def _vespagram_baz(stream, time_shift_table, starttime, endtime,
                       method="DLS", nthroot=1):
        """
        Estimating the azimuth or slowness vespagram.
        :type stream: stream: :class:`~obspy.core.stream.Stream`.
        :param stream: Stream object.
        :param time_shift_table: 2D timeshift table for each station in the
                                array. Each table gives the timeshift for all
                                slowness_x and slowness_y combinations.
        :type starttime: UTCDateTime
        :param starttime: Starttime of interest
        :type endtime: UTCDateTime
        :param endtime: Endtime of interest
        :type method: str
        :param method: Method used for computaion of the vespagram
                       Can be either:
                           -DlS: Delay and Sum
                           -PWS: Phase Weighted Stack
        :param nthroot: nth-root processing; nth gives the root (1,2,3,4),
         default 1 (no nth-root)
        :return: numpy.ndarray of beams with different slownesses
        """
        fs = stream[0].stats.sampling_rate
        if len(stream) != len(stream.select(sampling_rate=fs)):
            msg = 'All traces must have same sampling rate.'
            raise ValueError(msg)
        mini = min(value.min() for key, value in time_shift_table.items()
                   if key is not None)
        maxi = max(value.max() for key, value in time_shift_table.items()
                   if key is not None)
        spoint, _etime = _get_stream_offsets(stream, (starttime - mini),
                                             (endtime - maxi))

        # time shift table has slowness array under key `None`
        slownesses = time_shift_table[None]

        # Recalculate the maximum possible trace length
        ndat = int(((endtime - maxi) - (starttime - mini)) * fs)
        beams = np.zeros((len(slownesses), ndat), dtype='f8')

        stream.detrend()
        max_beam = 0.0
        slow = 0.0

        sll = slownesses[0]
        sls = slownesses[1] - sll

        # ids = [key for key in time_shift_table.keys() if key is not None]
        for _i, slowness in enumerate(slownesses):
            singlet = 0.0
            if method == 'DLS':
                for _j, tr in enumerate(stream.traces):
                    station = tr.id
                    s = spoint[_j] + int(
                        time_shift_table[station][_i] * fs + 0.5)
                    shifted = tr.data[s: s + ndat]
                    singlet += 1. / len(stream) * np.sum(shifted * shifted)
                    beams[_i] += 1. / len(stream) * \
                        np.power(np.abs(shifted), 1. / nthroot) * \
                        shifted / np.abs(shifted)

                beams[_i] = np.power(np.abs(beams[_i]), nthroot) * \
                    beams[_i] / np.abs(beams[_i])

                bs = np.sum(beams[_i] * beams[_i])
                bs /= singlet

                if bs > max_beam:
                    max_beam = bs
                    beam_max = _i
                    slow = slowness
                    if slow < 1e-8:
                        slow = 1e-8

            elif method == 'PWS':
                stack = np.zeros(ndat, dtype='c8')
                nstat = len(stream)
                for _j, tr in enumerate(stream.traces):
                    station = tr.id
                    s = spoint[_j] + int(
                        time_shift_table[station][_i] * fs + 0.5)
                    try:
                        shifted = sp.signal.hilbert(tr.data[s:s + ndat])
                    except IndexError:
                        break
                    phase = np.arctan2(shifted.imag, shifted.real)
                    stack.real += np.cos(phase)
                    stack.imag += np.sin(phase)
                coh = 1. / nstat * np.abs(stack)
                for _j, tr in enumerate(stream.traces):
                    station = tr.id
                    s = spoint[_j] + int(
                        time_shift_table[station][_i] * fs + 0.5)
                    shifted = tr.data[s: s + ndat]
                    singlet += 1. / nstat * np.sum(shifted * shifted)
                    beams[_i] += 1. / nstat * shifted * np.power(coh, nthroot)
                bs = np.sum(beams[_i] * beams[_i])
                bs = bs / singlet
                if bs > max_beam:
                    max_beam = bs
                    beam_max = _i
                    slow = np.abs(sll + _i * sls)
                    if slow < 1e-8:
                        slow = 1e-8
            else:
                msg = "Method '%s' unknown." % method
                raise ValueError(msg)

        return slow, beams, beam_max, max_beam

    @staticmethod
    def _geometry_dict_to_array(geometry):
        """
        Take a geometry dictionary (as provided by self.geometry, or by
        _get_geometry_xyz) and convert to a numpy array, as used in some
        methods.
        :type geometry: dict
        :param geometry: A dictionary with keys: SEED IDs and values:
                         dictionaries of 'latitude', 'longitude' and
                         'absolute_height_in_km'.
        :return: Sorted (by station) two dimensional numpy.ndarray of latitude,
                 longitude and height for each station.
        """
        geom_array = np.empty((len(geometry), 3))
        try:
            for _i, (key, value) in enumerate(sorted(list(geometry.items()))):
                geom_array[_i, 0] = value["x"]
                geom_array[_i, 1] = value["y"]
                geom_array[_i, 2] = value["z"]
        except KeyError:
            for _i, (key, value) in enumerate(sorted(list(geometry.items()))):
                geom_array[_i, 0] = float(value["longitude"])
                geom_array[_i, 1] = float(value["latitude"])
                geom_array[_i, 2] = value["absolute_height_in_km"]
        return geom_array

    def _correct_with_3dplane(self, geometry):
        """
        Correct a given array geometry with a best-fitting plane.

        :type geometry: dict
        :param geometry: A dictionary with keys: SEED IDs and values:
                         dictionaries of 'latitude', 'longitude' and
                         'absolute_height_in_km'.
        :return: The corrected geometry as dictionary, with the same keys as
            passed in.
        """
        # sort keys in the nested dict to be alphabetical:
        coord_sys_keys = sorted(list(geometry.items())[0][1].keys())
        if coord_sys_keys[0] == 'x':
            pass
        elif coord_sys_keys[0] == 'absolute_height_in_km':
            # set manually because order is important.
            coord_sys_keys = ['latitude', 'longitude', 'absolute_height_in_km']
        else:
            raise KeyError("Geometry dictionary does not have correct keys.")
        orig_geometry = geometry.copy()
        geo_a = self._geometry_dict_to_array(geometry)
        # compute barycenter of stations
        center = geo_a.sum(axis=0) / geo_a.shape[0]
        # compute basis and normal vector of the best fitting plane (vh[2])
        # the plane minimizes the squared distance of the points to the plane
        # taken from: "https://stackoverflow.com/questions/35070178/fit-plane-
        # to-a-set-of-points-in-3d-scipy-optimize-minimize-vs-scipy-linalg-lsts"
        u, s, vh = np.linalg.linalg.svd(geo_a-center)
        # satisfies the plane equation a*x + b*y + c*z = 0
        result = np.zeros((len(geometry), 3))
        # now we are seeking the station positions on that plane
        # geometry[:,2] += v[2,-1]
        n = vh[2, :]
        result[:, 0] = (geo_a[:, 0] - n[0] * (n[0] * geo_a[:, 0] +
                        geo_a[:, 1] * n[1] + n[2] * geo_a[:, 2]) /
                        (n[0] * n[0] + n[1] * n[1] + n[2] * n[2])**0.5)
        result[:, 1] = (geo_a[:, 1] - n[1] * (n[0] * geo_a[:, 0] +
                        geo_a[:, 1] * n[1] + n[2] * geo_a[:, 2]) /
                        (n[0] * n[0] + n[1] * n[1] + n[2] * n[2])**0.5)
        result[:, 2] = (geo_a[:, 2] - n[2] * (n[0] * geo_a[:, 0] +
                        geo_a[:, 1] * n[1] + n[2] * geo_a[:, 2]) /
                        (n[0] * n[0] + n[1] * n[1] + n[2] * n[2])**0.5) + \
                        center[2]
        geometry = result[:]
        print("Best fitting plane-coordinates :\n", geometry)

        # convert geometry array back to a dictionary.
        geodict = {}
        # The sorted list is necessary to match the station IDs (the keys in
        # the geometry dict) to the correct array row (or column?), same as is
        # done in _geometry_dict_to_array, but backwards.
        for _i, (key, value) in enumerate(sorted(
                list(orig_geometry.items()))):
            geodict[key] = {coord_sys_keys[0]: geometry[_i, 0],
                            coord_sys_keys[1]: geometry[_i, 1],
                            coord_sys_keys[2]: geometry[_i, 2]}
        geometry = geodict
        return geometry

    def show_distance_plot(self, stream, event, starttime, endtime,
                           plot_travel_times=True, vel_model='ak135'):
        """
        Plots distance dependent seismogramm sections.
        :param stream: Waveforms for the array processing.
        :type stream: :class:`obspy.core.stream.Stream`
        :param event: Earthquake position defined either by Obspy Event or
                      Obspy Origin class.
        :type event: :class:`~obspy.core.event.event.Event` or
            :class:`~obspy.core.event.origin.Origin`
        :param starttime: starttime of traces to be plotted
        :type starttime: UTCDateTime
        :param endtime: endtime of traces to be plotted
        :type endtime: UTCDateTime
        :param plot_travel_times: Flag weather phases are marked as traveltime
            plots in the section. obspy.taup is used to calculate the phases.
        :type: bool
        :param vel_model: 1D velocity model used for calculation of the
                          theoretical phase arrivals. A list of possible models
                          can be found here:
                          docs.obspy.org/packages/obspy.taup.html
        :type vel_model: str
        """
        model = TauPyModel(model=vel_model)
        stream = stream.slice(starttime=starttime, endtime=endtime).copy()

        if isinstance(event, Event):
            origin_ = event.origins[0]
        elif isinstance(event, Origin):
            origin_ = event
        else:
            raise TypeError("Only Obspy Events or Obspy Origins are supported"
                            " as event")

        event_depth_in_km = origin_.depth / 1000.0
        event_time = origin_.time

        self._attach_coords_to_stream(stream, origin_)

        cmap = plt.cm.get_cmap('jet')

        stream.traces = sorted(stream.traces,
                               key=lambda x: x.stats.coordinates.distance)[::-1]
        # One color for each trace.
        colors = [cmap(_i) for _i in np.linspace(0, 1, len(stream))]

        # Relative event times.
        times_array = stream[0].times() + \
            (stream[0].stats.starttime - event_time)

        distances = [tr.stats.coordinates.distance for tr in stream]
        min_distance = min(distances)
        max_distance = max(distances)
        distance_range = max_distance - min_distance
        stream_range = distance_range / 10.0

        # Normalize data and "shift to distance".
        stream.normalize()
        for tr in stream:
            tr.data *= stream_range
            tr.data += tr.stats.coordinates.distance

        plt.figure(figsize=(20, 12))
        for _i, tr in enumerate(stream):
            plt.plot(times_array, tr.data, label="%s.%s" % (tr.stats.network,
                                                            tr.stats.station),
                     color=colors[_i])
        plt.grid()
        plt.ylabel("Distance in degree to event")
        plt.xlabel("Time in seconds since event")
        plt.legend(loc="upper left")

        dist_min, dist_max = plt.ylim()

        if plot_travel_times:

            distances = defaultdict(list)
            ttimes = defaultdict(list)

            for i in np.linspace(dist_min, dist_max, 100):
                tts = model.get_travel_times(distance_in_degree=i,
                                             source_depth_in_km=event_depth_in_km)
                for arrival in tts:
                    name = arrival.name
                    distances[name].append(i)
                    ttimes[name].append(arrival.time)

            for key in distances.keys():
                min_distance = min(distances[key])
                max_distance = max(distances[key])
                min_tt_time = min(ttimes[key])
                max_tt_time = max(ttimes[key])

                if min_tt_time >= times_array[-1] or \
                        max_tt_time <= times_array[0] or \
                        (max_distance - min_distance) < \
                        0.8 * (dist_max - dist_min):
                    continue
                ttime = ttimes[key]
                dist = distances[key]
                if max(ttime) > times_array[0] + 0.9 * times_array.ptp():
                    continue
                plt.scatter(ttime, dist, s=0.5, zorder=-10, color="black",
                            alpha=0.8)
                plt.text(max(ttime) + 0.005 * times_array.ptp(),
                         dist_max - 0.02 * (dist_max - dist_min), key)

        plt.ylim(dist_min, dist_max)
        plt.xlim(times_array[0], times_array[-1])
        plt.title(event.short_str())

        plt.show()

    def align_phases(self, stream, event, phase_name, vel_model='ak135'):
        """
        Aligns all trace along a common theoretical phase arrival and returns
        the modified stream.

        :param stream: Stream containing all traces of the array.
        :type stream: :class:`~obspy.core.stream.Stream`.
        :param event: Event or Origin Object to calculate theoretical
                      travel-times from.
        :type event: :class:`~obspy.core.event.event.Event` or
                     :class:`~obspy.core.event.origin.Origin`
        :param phase_name: Phase by which the traces are aligned. Mostly
                           conventional naming, for more information see:
                            docs.obspy.org/packages/obspy.taup.html
        :type phase_name: str
        :param vel_model: 1D velocity model used for calculation of the
                          theoretical phase arrivals. A list of possible models
                          can be found here:
                          docs.obspy.org/packages/obspy.taup.html
        :type vel_model: str
        :return: Stream object with shifted traces.
        """
        starttime = max([tr.stats.starttime for tr in stream])
        stt = starttime
        endtime = min([tr.stats.endtime for tr in stream])
        stream.trim(starttime, endtime)
        shift = []
        stream = stream.copy()

        if isinstance(event, Event):
            origin_ = event.origins[0]
        elif isinstance(event, Origin):
            origin_ = event
        self._attach_coords_to_stream(stream, origin_)

        stream.traces = sorted(stream.traces,
                               key=lambda x: x.stats.coordinates.distance)[::-1]
        model = TauPyModel(model=vel_model)

        for tr in stream:
            tt = model.get_travel_times(
                distance_in_degree=tr.stats.coordinates.distance,
                source_depth_in_km=origin_.depth / 1000.0)
            if phase_name not in [arrival.name for arrival in tt]:
                raise ValueError('Selected Phase "{}" not occurring.'
                                 .format(phase_name))
            for t in tt:
                if t.name != phase_name:
                    continue
                tt_t = t.time
                break
            try:
                travel = origin_.time.timestamp + tt_t
                dtime = travel - stt.timestamp
                shift.append(dtime)
            except:
                break
        shift = np.asarray(shift)
        shift -= shift[0]
        self.shifttrace_freq(stream, -shift)
        return stream

    @staticmethod
    def shifttrace_freq(stream, t_shift):
        if isinstance(stream, Stream):
            for i, tr in enumerate(stream):
                ndat = tr.stats.npts
                samp = tr.stats.sampling_rate
                nfft = next_pow_2(ndat)
                nfft *= 2
                tr1 = np.fft.rfft(tr.data, int(nfft))
                for k in range(0, int(nfft / 2), 1):
                    tr1[k] *= complex( np.cos((t_shift[i] * samp) * (k / float(nfft)) * 2. * np.pi) \
                               -np.sin((t_shift[i] * samp) * (k / float(nfft)) * 2. * np.pi))

                tr1 = np.fft.irfft(tr1, nfft)
                tr.data = tr1[0:ndat]


if __name__ == '__main__':
    import doctest

    doctest.testmod(exclude_empty=True)
