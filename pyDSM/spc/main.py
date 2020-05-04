import numpy as np
import spctime
import sys
try:
    import tish
except ModuleNotFoundError:
    sys.path.append('../../lib')
    import tish

import matplotlib.pyplot as plt

def plot(u, sampling_hz):
    nr = u.shape[1]
    npts = u.shape[2]
    fig, axes = plt.subplots(nr, 1)
    ts = np.array(list(range(npts))) / sampling_hz
    for ir, ax in enumerate(axes.ravel()):
        ax.plot(ts, u[2, ir])
        ax.set(xlabel='Time (s)',
            ylabel='Ground vel.')

if __name__ == '__main__':
    inputs = tish.pinput_fromfile('../../src_f90/tish/example/dsm_accuracy_check/AK135_SH.inf')
    spcs = tish.tish(*inputs, False)

    tlen = inputs[3]
    nspc = inputs[4]
    omegai = inputs[5]
    sampling_hz = 20

    spct = spctime.SpcTime(tlen, nspc, sampling_hz, omegai)

    u = spct.spctime(spcs)

    plot(u, sampling_hz)
    plt.show()
