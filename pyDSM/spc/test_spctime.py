import spctime
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

try:
    import tish
except ModuleNotFoundError:
    sys.path.append('../../lib')
    import tish

rootdsm = os.path.abspath('../../src_f90/tish/example/dsm_accuracy_check/')

def get_u_pydsm():
    inputs = tish.pinput_fromfile(os.path.join(rootdsm, 'AK135_SH.inf'))
    spcs = tish.tish(*inputs, False)

    tlen = inputs[3]
    nspc = inputs[4]
    omegai = inputs[5]
    sampling_hz = 20

    spct = spctime.SpcTime(tlen, nspc, sampling_hz, omegai)

    u = spct.spctime(spcs)

    return u, tlen

def get_u_dsm():
    u = np.loadtxt(os.path.join(rootdsm, 'sac/109C_TA.200702131456A.T.txt'))
    return u

def plot(udsm, upydsm, tlen):
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(10,8))
    ts = np.linspace(0, tlen, len(udsm))
    ax0.plot(ts, udsm, label='dsm')
    ax0.plot(ts, upydsm, label='pydsm', color='red')

    residuals = upydsm - udsm
    ax1.plot(ts, residuals, color='blue', label='residuals')

    ax0.legend()
    ax0.set(xlabel='Time (s)', ylabel='Ground velocity', xlim=[0,3276.8])
    ax1.legend()
    ax1.set(xlabel='Time (s)', ylabel='Residuals', xlim=[0,3276.8])

    return fig, (ax0, ax1)

if __name__ == '__main__':
    udsm = get_u_dsm()

    print('Computing waveform using pyDSM (np=512)')
    upydsms, tlen = get_u_pydsm()
    upydsm = upydsms[2, 0]

    filename = './waveform_accuracy.pdf'
    print('Saving DSM and pyDSM waveform comparison in {}'.format(filename))
    _ = plot(udsm, upydsm, tlen)
    plt.savefig(filename, bbox_inches='tight')

    assert(np.allclose(udsm, upydsm, rtol=1e-10))
    print('All passed!')