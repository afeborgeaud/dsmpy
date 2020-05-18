import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pydsm import dsm, rootdsm_sh
from pydsm.spc import spctime
from pydsm.seismicmodel import SeismicModel

from pydsm._tish import _pinput, _tish

def get_input():
    seismic_model = SeismicModel.ak135()

    mt = dsm.MomentTensor(Mrr=-0.193, Mrt=0.0998, Mrp=-0.148,
                                 Mtt=0.0571, Mtp=-0.141, Mpp=0.136)
    source_time_function = None

    event = dsm.Event('', -1.600000023841858, -78.05000305175781,
                      6371-6195.2, mt.to_array(),
                      source_time_function=source_time_function)
    stations = [
        dsm.Station('109C', 'TA', 32.88890075683594, -117.1051025390625)]
    tlen = 3276.8
    nspc = 64
    source_time_function = None
    sampling_hz = 20

    pydsm_input = dsm.PyDSMInput.input_from_arrays(
        event, stations, seismic_model,
        tlen, nspc, sampling_hz)
    return pydsm_input

def get_u_pydsm():
    pydsm_input = get_input()

    outputs = dsm.compute(pydsm_input)
    outputs.to_time_domain()
    return outputs.us, outputs.ts

def get_u_dsm():
    u = np.loadtxt(os.path.join(rootdsm_sh,
                                'sac_64/109C_TA.200702131456A.T.txt'))
    return u


def plot(ts, udsm, upydsm):
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 8))
    ax0.plot(ts, udsm, label='dsm')
    ax0.plot(ts, upydsm, label='pydsm', color='red')

    residuals = upydsm - udsm
    ax1.plot(ts, residuals, color='blue', label='residuals')

    ax0.legend()
    ax0.set(xlabel='Time (s)', ylabel='Ground velocity',
            xlim=[0, 3276.8])
    ax1.legend()
    ax1.set(xlabel='Time (s)', ylabel='Residuals', xlim=[0, 3276.8])

    return fig, (ax0, ax1)


if __name__ == '__main__':
    udsm = get_u_dsm()

    print('Computing waveform using pyDSM')
    upydsms, ts = get_u_pydsm()
    upydsm = upydsms[2, 0]

    filename = 'figures/waveform_accuracy.pdf'
    print('Saving DSM and pyDSM waveform comparison in {}'
          .format(filename))
    _ = plot(ts, udsm, upydsm)
    plt.savefig(filename, bbox_inches='tight')

    assert (np.allclose(udsm, upydsm, rtol=1e-10))
    print('All passed!')
