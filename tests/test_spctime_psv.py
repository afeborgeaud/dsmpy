import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pydsm import dsm, rootdsm_psv
from pydsm.spc import spctime

from pydsm._tish import _pinput, _tish


def get_u_pydsm():
    parameter_file = os.path.join(rootdsm_psv, 'test1.inf')
    inputs = dsm.PyDSMInput.input_from_file(
        parameter_file, sampling_hz=20,
        source_time_function=None, mode=0)
    outputs_sh, outputs_psv = dsm.compute(inputs)
    outputs_sh.to_time_domain()
    outputs_psv.to_time_domain()
    return outputs_sh.us, outputs_sh.ts, outputs_psv.us


def get_u_dsm():
    u = np.loadtxt(os.path.join(rootdsm_psv,
                                'sac_64/109C_TA.200702131456A.T.txt'))
    return u


def plot(ts, udsm, upydsm):
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 8))
    ax0.plot(ts, udsm, label='dsm')
    #ax0.plot(ts, upydsm, label='pydsm', color='red')

    residuals = upydsm #- udsm
    ax1.plot(ts, residuals, color='blue', label='residuals')

    ax0.legend()
    ax0.set(xlabel='Time (s)', ylabel='Ground velocity')
            #xlim=[0, 3276.8])
    ax1.legend()
    ax1.set(xlabel='Time (s)', ylabel='Residuals')
        #, xlim=[0, 3276.8])

    return fig, (ax0, ax1)


if __name__ == '__main__':
    #udsm = get_u_dsm() # TODO add psv file

    print('Computing waveform using pyDSM')
    upydsms, ts, upydsm_psv = get_u_pydsm()
    upydsm = upydsms[1, 0]
    upydsm_psv = upydsm_psv[1, 0]

    filename = 'figures/waveform_accuracy_psv.pdf'
    print('Saving DSM and pyDSM waveform comparison in {}'
          .format(filename))
    #_ = plot(ts, udsm, upydsm) #TODO add psv file
    _ = plot(ts, upydsm, upydsm_psv)
    plt.savefig(filename, bbox_inches='tight')

    #assert (np.allclose(udsm, upydsm, rtol=1e-10))
    #print('All passed!')
