"""Test dsmpy comparing time-domain waveforms computed using dsmpy
with spctime (freqency to time domain)
to waveforms computed using DSM and the Java tool Kibrary
(https://github.com/kensuke1984/Kibrary)"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from dsmpy import dsm, rootdsm_sh
from dsmpy.spc import spctime
from dsmpy._tish import _pinput, _tish


def _get_u_pydsm():
    parameter_file = os.path.join(rootdsm_sh, 'AK135_SH_64.inf')
    inputs = dsm.PyDSMInput.input_from_file(
        parameter_file, sampling_hz=20,
        source_time_function=None, file_mode=2)
    outputs = dsm.compute(inputs, mode=2)
    outputs.to_time_domain()
    return outputs


def _get_u_dsm():
    u = np.loadtxt(os.path.join(rootdsm_sh,
                                'sac_64/109C_TA.200702131456A.T.txt'))
    return u


def _plot(ts, udsm, upydsm):
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


def test_spctime_sh():
    udsm = _get_u_dsm()

    outputs = _get_u_pydsm()
    upydsms, ts = outputs.us, outputs.ts
    upydsm = upydsms[2, 0]
    
    print(udsm[:5])
    print(upydsm[:5])

    assert (np.allclose(udsm, upydsm, rtol=1e-10))

    return outputs, udsm


if __name__ == '__main__':
    outputs, udsm = test_spctime_sh()

    # Plot DSM and pydsm waveform comparison
    if 'figures' not in os.listdir('.'):
        os.mkdir('./figures')
    filename = 'figures/waveform_accuracy_sh.pdf'
    _plot(outputs.ts, udsm, outputs.us[2, 0])
    plt.savefig(filename, bbox_inches='tight')

    # Write to SAC files
    root_path = 'figures'
    outputs.write(root_path, 'sac')
