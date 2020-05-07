import numpy as np
from pydsm import dsm
from pydsm.spc.spctime import SourceTimeFunction
from pydsm import rootdsm
import matplotlib.pyplot as plt
from scipy.integrate import trapz


if __name__ == '__main__':
    sampling_hzs = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20]
    #samplingHzs = [0.5, 1, 2, 5, 10, 20, 50, 100]
    parameter_file = rootdsm + '/AK135_SH_64.inf'

    Es = []
    for sampling_hz in sampling_hzs:
        inputs = dsm.PyDSMInput.input_from_file(parameter_file, 
            sampling_hz=sampling_hz)
        outputs = dsm.compute(inputs)
        outputs.to_time_domain()
        u = outputs.u[2,0]
        E = trapz(u**2, dx=1/sampling_hz**2)
        Es.append(E)
    
    fig, ax = plt.subplots(1)
    ax.plot(sampling_hzs, Es, '-o', color='red')
    ax.set(xlabel='Sampling freq. (Hz)',
        ylabel='Energy')
    plt.savefig('figures/test_samplingHz_64.pdf', bbox_inches='tight')