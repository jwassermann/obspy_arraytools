import glob
import obspy
import os

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

folders = ["01_GRSN_GRF", "02_MM", "03_YKA", "04_KNET", "05_EAGLE"]


for path in folders:
    print ""
    print 80 * "="
    print 80 * "="
    print "Testing folder '%s'..." % path

    output_path = os.path.join(path, "_test_output")

    if os.path.exists(output_path):
        msg = ("Folder '%s' has already been tested. Remove '_test_output' "
               "to test again." % path)
        print msg
        continue
    os.makedirs(output_path)

    mseed = glob.glob(os.path.join(path, "*.mseed"))
    station_xml = glob.glob(os.path.join(path, "*.xml"))
    quakeml = glob.glob(os.path.join(path, "*.qml"))

    if len(mseed) != 1 and "02" not in path:
        msg = "ERROR in folder '%s': Not exactly one MiniSEED file!" % path
        print msg
        print "Will be skipped!"
        continue
    if len(station_xml) != 1:
        msg = "ERROR in folder '%s': Not exactly one StationXML file!" % path
        print msg
        print "Will be skipped!"
        continue

    if len(quakeml) != 1:
        msg = "ERROR in folder '%s': Not exactly one QuakeML file!" % path
        print msg
        print "Will be skipped!"
        continue

        raise Exception(msg)

    if "02" not in path:
        mseed = mseed[0]
        station_xml = station_xml[0]
        quakeml = quakeml[0]

        st = obspy.read(mseed)
        inv = obspy.read_inventory(station_xml, format="stationxml")
        cat = obspy.readEvents(quakeml)

        ev = cat[0]
        org = ev.preferred_origin() or ev.origins[0]
        ev_lat = org.latitude
        ev_lon = org.longitude

        st.attach_response(inv)

        for tr in st:
            if not hasattr(tr.stats, "response") or not tr.stats.response:
                msg = "WARNING: No response for channel '%s'." % tr.id
                print msg

        st.plot(outfile=os.path.join(output_path, "uncorrected.png"),
                equal_scale=False)

        st.remove_response(output="VEL", pre_filt=(0.01, 0.1, 10, 20))
        st.plot(outfile=os.path.join(output_path, "corrected.png"))
        required_stations = set([".".join(tr.id.split(".")[:2]) for tr in st])

    else:
        mseed = glob.glob(os.path.join(path, "*.227"))
        station_xml = station_xml[0]
        quakeml = quakeml[0]
        inv = obspy.read_inventory(station_xml, format="stationxml")
        cat = obspy.readEvents(quakeml)
        ev = cat[0]
        org = ev.preferred_origin() or ev.origins[0]
        ev_lat = org.latitude
        ev_lon = org.longitude

        required_stations = set()

        for filename in mseed:
            print filename
            st = obspy.read(filename)
            if len(st) != 1:
                msg = "File '%s' contains %i traces!!!" % \
                    (filename, len(st))
                print msg
            required_stations.add(".".join(st[0].id.split(".")[:2]))

            st.attach_response(inv)

            for tr in st:
                if not hasattr(tr.stats, "response") or not tr.stats.response:
                    msg = "WARNING: No response for channel '%s'." % tr.id
                    print msg

            st.plot(outfile=os.path.join(output_path,
                    "%s_uncorrected.png" % os.path.basename(filename)),
                    equal_scale=False)

            st.remove_response(output="VEL", pre_filt=(0.01, 0.1, 10, 20))
            st.plot(outfile=os.path.join(
                output_path, "%s_corrected.png" % os.path.basename(filename)))

    coords = {}

    for network in inv:
        network_id = network.code
        for station in network:
            station_id = station.code
            id = "%s.%s" % (network_id, station_id)
            if id not in required_stations:
                continue
            lat = station.latitude
            lon = station.longitude

            if lat is None or lon is None:
                continue

            coords[id] = {"latitude": float(station.latitude),
                          "longitude": float(station.longitude)}

    available_stations = set(coords.keys())

    if not required_stations.issubset(available_stations):
        missing = required_stations.difference(available_stations)
        msg = "WARNINGS: Coordinates missing for stations %s" % missing
        print msg
        continue

    lons = []
    lats = []
    for key in required_stations:
        a = coords[key]
        lons.append(a["longitude"])
        lats.append(a["latitude"])

    plt.figure(figsize=(30, 20))
    m = Basemap(projection='hammer', lon_0=180)
    x, y = m(lons, lats)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966', lake_color='#99ffff')
    m.scatter(x, y, 30, marker='o', color='k', zorder=100000)
    plt.title("Station locations")

    for lon, lat in zip(lons, lats):
        m.drawgreatcircle(ev_lon, ev_lat, lon, lat, lw=2, zorder=10000000)

    x, y = m((ev_lon,), (ev_lat,))
    m.scatter(x, y, 50, marker="o", color="red", zorder=100000000000)

    plt.savefig(os.path.join(output_path, "global_map.png"))

    min_lon = min(lons)
    max_lon = max(lons)
    min_lat = min(lats)
    max_lat = max(lats)

    if "03" in path:
        min_lat -= 1.0
        max_lat += 1.0
        min_lon -= 1.0
        max_lon += 1.0
    else:
        min_lat -= 10.0
        max_lat += 10.0
        min_lon -= 10.0
        max_lon += 10.0

    plt.figure(figsize=(30, 20))
    m = Basemap(projection='cyl', llcrnrlat=min_lat, urcrnrlat=max_lat,
                llcrnrlon=min_lon, urcrnrlon=max_lon, resolution='i')
    x, y = m(lons, lats)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966', lake_color='#99ffff')
    m.drawcountries(linewidth=0.2, zorder=9000)
    m.scatter(x, y, 10, marker='o', color='k', zorder=100000)
    plt.title("Station locations")

    for lon, lat in zip(lons, lats):
        m.drawgreatcircle(ev_lon, ev_lat, lon, lat, lw=2, zorder=10000000)

    x, y = m((ev_lon,), (ev_lat,))
    m.scatter(x, y, 50, marker="o", color="red", zorder=100000000000)

    plt.savefig(os.path.join(output_path, "local_map.png"))
