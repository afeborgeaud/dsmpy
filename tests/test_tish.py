import numpy as np
from dsmpy._tish import _tish, _pinput
from dsmpy import rootdsm_sh
import time
import os


def test_tish():
    # parameter_file = rootdsm_sh + 'AK135_SH_64.inf'
    parameter_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'input_files/AK135_SH_64.inf')
    inputs = _pinput(parameter_file)
    write_to_file = False
    start = time.time()
    u = _tish(*inputs, write_to_file)
    u2 = _tish(*inputs, write_to_file)
    end = time.time()
    print('tish finished in {} s'.format(end - start))
    error_re = abs(u[2, 0, -1].real + 4.7034875e-10) / 4.7034875e-10
    error_im = abs(u[2, 0, -1].imag + 1.7664374e-11) / 1.7664374e-11

    # until now, 1e-7 works for -O3 optimisation
    # using -fast-math gfortran flag (gcc 4.4.7)
    # gives and error_im of ~1.4e-7, so we changed
    # to 1e-6, which is still a relatively small error
    assert error_re < 1e-6
    assert error_im < 1e-6
    # to test value leaks when calling _tish 2 times
    assert np.allclose(u, u2)
    print("All passed")


if __name__ == '__main__':
    test_tish()
