import numpy as np
from pyDSM import dsm
from pyDSM import rootdsm
from pyDSM.spc.spctime import SourceTimeFunction

import matplotlib.pyplot as plt

def plot(outputs):
    nr = outputs.get_nr()
    fig, axes = plt.subplots(nr, 1, sharex=True)
    for ir, ax in enumerate(axes.ravel()):
        ax.plot(outputs.ts, outputs.u[2, ir], label=outputs.stations[ir])
        ax.set(ylabel='Ground vel.')
        ax.legend()
    axes[-1].set_xlabel('Time (s)')
    fig.suptitle(outputs.event.eventID)
    return fig, axes

if __name__ == '__main__':
    parameter_file = rootdsm + '/AK135_SH_64.inf'
    half_duration = 2.

    inputs = dsm.PyDSMinput(parameter_file)
    stf = SourceTimeFunction.triangle(half_duration, inputs)
    inputs.set_sourcetimefunction(stf)
    outputs_4 = dsm.compute(inputs)
    outputs_4.to_time_domain()

    inputs = dsm.PyDSMinput(parameter_file, samplingHz=20)
    stf = SourceTimeFunction.triangle(half_duration, inputs)
    inputs.set_sourcetimefunction(stf)
    outputs_20 = dsm.compute(inputs)
    outputs_20.to_time_domain()

    print(len(outputs_4.ts), len(outputs_20.ts))

    fig, axes = plot(outputs_4)
    for ir, ax in enumerate(axes.ravel()):
        ax.plot(outputs_20.ts, outputs_20.u[2, ir], color='red')
    plt.show()
