import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from dsmpy import dsm, rootdsm_psv
from dsmpy.spc import spctime
from dsmpy._tish import _pinput, _tish


def _get_u_pydsm():
    parameter_file = os.path.join(rootdsm_psv, 'test11.inf')
    inputs = dsm.PyDSMInput.input_from_file(
        parameter_file, sampling_hz=20,
        source_time_function=None, file_mode=1)
    outputs = dsm.compute(inputs, mode=1)
    outputs.to_time_domain()
    return outputs


def _get_u_dsm():
    u = np.loadtxt(
        os.path.join(
            rootdsm_psv, 'sac_64/ABC_DSM.200702131456A.Zs.txt'))
    return u[:, 1]


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


def test_spctime_psv():
    udsm = _get_u_dsm()

    outputs = _get_u_pydsm()
    upydsms, ts = outputs.us, outputs.ts
    upydsm = upydsms[0, 0]  # Z-component, first station

    assert (np.allclose(udsm, upydsm, rtol=1e-10))

    return outputs, udsm


if __name__ == '__main__':
    outputs, udsm = test_spctime_psv()

    if 'figures' not in os.listdir('.'):
        os.mkdir('./figures')
    filename = 'figures/waveform_accuracy_psv.pdf'
    _plot(outputs.ts, udsm, outputs.us[0, 0])
    plt.savefig(filename, bbox_inches='tight')

# dir = '/Users/navy/git/dsmpy/dsmpy/src_f90/tipsv/examples/sac_64/200702131456A'
# paths = os.listdir(dir)
# paths = [os.path.join(dir, x) for x in paths]
# traces = [read(p) for p in paths]
# for i, trace in enumerate(traces):
#     tr = trace[0]
#     dt = tr.stats.delta
#     us = tr.data
#     ts = np.array([i * dt for i in range(len(us))])
#     outpath = paths[i] + '.txt'
#     with open(outpath, 'w') as f:
#         for t, u in zip(ts, us):
#             f.write(f'{t} {u}\n')
