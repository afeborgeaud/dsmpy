from pydsm import rootdsm_psv, rootdsm_sh
from pydsm import dsm
import os
import matplotlib.pyplot as plt
from pydsm.windows import Windows

if __name__ == '__main__':
    # compute P-SV + SH using P-SV input file
    parameter_file_psv = os.path.join(rootdsm_psv, 'test4.inf')
    input = dsm.PyDSMInput.input_from_file(
        parameter_file_psv, sampling_hz=20,
        source_time_function=None, mode=1)
    output = dsm.compute(input, mode=0)
    
    windows = Windows(
        output.event, output.stations, 'prem', ['ScS', 'Sdiff'])
    window_width = 80.
    output_windowed = output.window_spcs(windows, window_width)

    output.to_time_domain()
    output_windowed.to_time_domain()

    fig, axes = output.plot()
    _, axes = output_windowed.plot(axes=axes)

    fig_spc, axes_spc = output.plot_spcs()
    _, axes_spc = output_windowed.plot(axes=axes_spc)
    
    plt.show()