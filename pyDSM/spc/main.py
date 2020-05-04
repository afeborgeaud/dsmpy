import numpy as np
from pyDSM import dsm
from pyDSM import rootdsm
from pyDSM.spc.spctime import SourceTimeFunction

import matplotlib.pyplot as plt

def plot(outputs):
    nr = outputs.get_nr()
    fig, axes = plt.subplots(nr, 1, sharex=True)
    for ir, ax in enumerate(axes.ravel()):
        ax.plot(ts, outputs.u[2, ir], label=outputs.stations[ir])
        ax.set(ylabel='Ground vel.')
        ax.legend()
    axes[-1].set_xlabel('Time (s)')
    fig.suptitle(outputs.event.eventID)

if __name__ == '__main__':
    parameter_file = rootdsm + '/AK135_SH.inf'
    half_duration = 10.
    inputs = dsm.DSMinput(parameter_file)
    stf = SourceTimeFunction.triangle(half_duration, inputs)
    inputs.set_sourcetimefunction(stf)
    outputs = dsm.compute(inputs)
    outputs.to_time_domain()

    u = outputs.u
    ts = outputs.ts

    plot(outputs)
    plt.show()
