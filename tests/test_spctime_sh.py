import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pydsm import dsm, rootdsm_sh
from pydsm.spc import spctime
from pydsm._tish import _pinput, _tish


def get_u_pydsm():
    parameter_file = os.path.join(rootdsm_sh, 'AK135_SH_64.inf')
    inputs = dsm.PyDSMInput.input_from_file(
        parameter_file, sampling_hz=20,
        source_time_function=None, mode=2)
    outputs = dsm.compute(inputs, mode=2)
    outputs.to_time_domain()
    return outputs


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
    outputs = get_u_pydsm()
    upydsms, ts = outputs.us, outputs.ts
    upydsm = upydsms[2, 0]

    filename = 'figures/waveform_accuracy_sh.pdf'
    print('Saving DSM and pyDSM waveform comparison in {}'
          .format(filename))
    _ = plot(ts, udsm, upydsm)
    plt.savefig(filename, bbox_inches='tight')

    assert (np.allclose(udsm, upydsm, rtol=1e-10))
    print('All passed!')

    print('Write to SAC')
    root_path = 'figures'
    outputs.write(root_path, 'sac')
