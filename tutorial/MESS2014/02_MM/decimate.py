import glob
import numpy as np
import obspy
import os
from scipy import signal


def zerophase_chebychev_lowpass_filter(trace, freqmax):
    """
    Custom Chebychev type two zerophase lowpass filter useful for
    decimation filtering.

    This filter is stable up to a reduction in frequency with a factor of
    10. If more reduction is desired, simply decimate in steps.

    Partly based on a filter in ObsPy.

    :param trace: The trace to be filtered.
    :param freqmax: The desired lowpass frequency.
    """
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (trace.stats.sampling_rate * 0.5)  # stop band frequency
    wp = ws  # pass band frequency

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = signal.cheby2(order, rs, wn, btype="low", analog=0, output="ba")

    # Apply twice to get rid of the phase distortion.
    trace.data = signal.filtfilt(b, a, trace.data)

FOLDER = os.path.join("other", "full_sampling_rate")

for filename in glob.glob(os.path.join(FOLDER, "*.227")):
    print "Decimating file '%s' ..." % filename
    st = obspy.read(filename)

    for tr in st:
        decimation_factor = 5
        new_nyquist = tr.stats.sampling_rate / 2.0 / float(decimation_factor)
        zerophase_chebychev_lowpass_filter(tr, new_nyquist)
        tr.decimate(factor=decimation_factor, no_filter=True)
        tr.data = np.require(tr.data, "float32")


    st.write(os.path.basename(filename), format="mseed")
