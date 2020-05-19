import numpy as np
from numpy import f2py
import sys
import subprocess
import os
import glob
import shutil
from pydsm._tipsv import _tipsv, _pinput
from pydsm import rootdsm_psv
import time


def test_tipsv():
    parameter_file = rootdsm_psv + 'test1.inf'
    inputs = _pinput(parameter_file)
    write_to_file = False
    start = time.time()
    u = _tipsv(*inputs, write_to_file)
    u2 = _tipsv(*inputs, write_to_file)
    end = time.time()
    print('tipsv finished in {} s'.format(end - start))
    error_re = abs(
        u[2, 0, -1].real + 1.54309401723643564E-011) /-1.54309401723643564E-011
    error_im = abs(
        u[2, 0, -1].imag - 5.97994400276396742E-012) / 5.97994400276396742E-012
    assert error_re < 1e-7
    assert error_im < 1e-7
    # to test value leaks when calling _tipsv 2 times
    assert np.allclose(u, u2)
    print("All passed")


if __name__ == '__main__':
    test_tipsv()
