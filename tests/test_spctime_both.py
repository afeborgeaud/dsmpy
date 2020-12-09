from dsmpy import rootdsm_psv, rootdsm_sh
from dsmpy import dsm
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # compute P-SV using P-SV input file
    parameter_file_psv = os.path.join(rootdsm_psv, 'test4.inf')
    inputs1 = dsm.PyDSMInput.input_from_file(
        parameter_file_psv, sampling_hz=20,
        source_time_function=None, mode=1)
    outputs_psv = dsm.compute(inputs1, mode=1)
    outputs_psv.to_time_domain()

    # compute SH using SH input file
    parameter_file_sh = os.path.join(rootdsm_sh, '../test4.inf')
    inputs2 = dsm.PyDSMInput.input_from_file(
        parameter_file_sh, sampling_hz=20,
        source_time_function=None, mode=2)
    outputs_sh = dsm.compute(inputs2, mode=2)
    outputs_sh.to_time_domain()

    # compute P-SV + SH using P-SV input file
    parameter_file_psv = os.path.join(rootdsm_psv, 'test4.inf')
    inputs3 = dsm.PyDSMInput.input_from_file(
        parameter_file_psv, sampling_hz=20,
        source_time_function=None, mode=1)
    outputs_psv_sh = dsm.compute(inputs3, mode=0)
    outputs_psv_sh.to_time_domain()

    # plot
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
    ax1.plot(outputs_sh.ts, outputs_sh.us[2, 0], label='T SH',
        color='blue')
    ax1.plot(
        outputs_sh.ts, outputs_psv.us[2, 0],
        label='T P-SV', color='red')
    ax2.plot(
        outputs_sh.ts, outputs_sh.us[1, 0], label='R SH',
        color='blue')
    ax2.plot(
        outputs_sh.ts, outputs_psv.us[1, 0],
        label='R P-SV', color='red')
    ax3.plot(outputs_sh.ts, outputs_psv_sh.us[2, 0], label='T P-SV+SH',
        color='green')
    ax4.plot(
        outputs_sh.ts, outputs_psv_sh.us[1, 0],
        label='R P-SV+SH', color='green')
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax4.set(xlabel='Time (s)', ylabel='Velocity ()')
    
    plt.savefig('figures/waveform_accuracy_both.pdf')