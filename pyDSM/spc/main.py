import numpy as np
from pyDSM import dsm
from pyDSM import rootdsm

import matplotlib.pyplot as plt

def plot(outputs):
    nr = outputs.get_nr()
    fig, axes = plt.subplots(nr, 1)
    for ir, ax in enumerate(axes.ravel()):
        ax.plot(ts, u[2, ir])
        ax.set(xlabel='Time (s)',
            ylabel='Ground vel.')

if __name__ == '__main__':
    parameter_file = rootdsm + '/AK135_SH_64.inf'
    inputs = dsm.DSMinput(parameter_file)
    outputs = dsm.compute(inputs)
    outputs.to_time_domain()

    u = outputs.u
    ts = outputs.ts

    plot(outputs)
    plt.show()
