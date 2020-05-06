import numpy as np
from pyDSM import dsm
from pyDSM.spc.spctime import SourceTimeFunction
from pyDSM import rootdsm
import matplotlib.pyplot as plt
from scipy.integrate import trapz


if __name__ == '__main__':
    samplingHzs = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20]
    #samplingHzs = [0.5, 1, 2, 5, 10, 20, 50, 100]
    parameter_file = rootdsm + '/AK135_SH_64.inf'

    Es = []
    for samplingHz in samplingHzs:
        inputs = dsm.DSMinput(parameter_file, samplingHz=samplingHz)
        outputs = dsm.compute(inputs)
        outputs.to_time_domain()
        u = outputs.u[2,0]
        E = trapz(u**2, dx=1/samplingHz**2)
        Es.append(E)
    
    fig, ax = plt.subplots(1)
    ax.plot(samplingHzs, Es, '-o', color='red')
    ax.set(xlabel='Sampling freq. (Hz)',
        ylabel='Energy')
    plt.savefig('figures/test_samplingHz_64.pdf', bbox_inches='tight')