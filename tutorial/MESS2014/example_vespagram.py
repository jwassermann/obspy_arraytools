#!/usr/bin/env python

import os
from glob import glob
from obspy import read, UTCDateTime, read_inventory, readEvents
from mess2014 import vespagram

# Morocco Array
stt = UTCDateTime("2012-08-14T03:15:30.0")
e = UTCDateTime("2012-08-14T03:25:00.0")
folder = '02_MM/'
active_channel = 'HHZ'
# GRSN Array
#stt = UTCDateTime("1991-12-17T06:49:50.0")
#e = UTCDateTime("1991-12-17T06:50:30.0")
#folder = '01_GRSN_GRF/'
#active_channel = 'BHZ'
# Yellowknife Array
#stt = UTCDateTime("2012-08-14T03:07:30.0")
#e = UTCDateTime("2012-08-14T03:10:00.0")
#folder = '03_YKA/'
#active_channel = 'SHZ'
# KNET Array
#stt = UTCDateTime("2012-08-14T03:07:30.0")
#e = UTCDateTime("2012-08-14T03:10:00.0")
#folder = '04_KNET/'
#active_channel = 'BHZ'

# frequency range for array analysis
frqlow = 0.1
frqhigh = 1.0
# filter waveforms?
filter = True
# static topography correction?
static3D = False
# correction velocity for static topography correction
vc = 4.8
# method to use for array analysis ("FK", "DLS")
method = 'DLS'

# slowness min/max/step in s/deg
sl = (0, 10, 0.1)

plot_trace = True
baz = 17.
scale = 10
nthroot = 4
align = False


mseedfile = glob(os.path.join(folder, "??.mseed"))
stationxmlfile = glob(os.path.join(folder, "??.xml"))
quakemlfile = glob(os.path.join(folder, "??.qml"))

stream = read(mseedfile[0])
stream = stream.select(channel=active_channel)
stream.trim(stt, e)

inventory = read_inventory(stationxmlfile[0], format="STATIONXML")
event = readEvents(quakemlfile[0])[0]

vespagram(stream, event, inventory, method=method, frqlow=frqlow,
          frqhigh=frqhigh, baz=baz, scale=scale, nthroot=nthroot,
          filter=filter, static3D=static3D, vel_corr=vc, sl=sl)
