import spctime
import numpy as np
import tish
import matplotlib.pyplot as plt

def get_u_pydsm():
    inputs = tish.pinput_fromfile('../example/AK135_SH.inf') #'../example/dsm_accuracy_check/AK135_SH.inf'
    spc = tish.tish(*inputs, False)

    tlen = inputs[3]
    nspc = inputs[4]
    omegai = inputs[5]
    sampling_hz = 20

    spct = spctime.SpcTime(tlen, nspc, sampling_hz, omegai)

    spc3 = spc[2, 0 , :]
    u = spct.spctime(spc3)

    return u, tlen

def get_u_dsm():
    u = np.loadtxt('../example/dsm_accuracy_check/sac/109C_TA.200702131456A.T.txt')
    return u

def plot(udsm, upydsm, tlen):
    fig, (ax0, ax1) = plt.subplots(2,1)
    ts = np.linspace(0, tlen, len(udsm))
    ax0.plot(ts, udsm, label='dsm')
    ax1.plot(ts, upydsm, label='pydsm')
    ax0.legend()
    ax1.legend()
    return fig, (ax0, ax1)

if __name__ == '__main__':
    udsm = get_u_dsm()
    upydsm, tlen = get_u_pydsm()

    _ = plot(udsm, upydsm, tlen)
    plt.show()