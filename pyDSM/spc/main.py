import numpy as np
import spctime
import tish

import matplotlib.pyplot as plt

def plot(u, sampling_hz):
    fig, ax = plt.subplots(1, 1)
    ts = np.array(list(range(len(u)))) / sampling_hz
    ax.plot(ts, u)
    ax.set(xlabel='Time (s)',
        ylabel='Ground velocity')
    #plt.show()

if __name__ == '__main__':
    inputs = tish.pinput_fromfile('../AK135_SH.inf')
    spc = tish.tish(*inputs, False)

    tlen = inputs[3]
    nspc = inputs[4]
    omegai = inputs[5]
    sampling_hz = 20
    print(omegai, tlen, nspc)

    spctime = spctime.SpcTime(tlen, nspc, sampling_hz, omegai)

    spc3 = spc[2, 0 , :]

    u = spctime.spctime(spc3)

    plt.plot(np.abs(spc3))
    plot(u, sampling_hz)

    plt.show()
