import numpy as np
from numpy import f2py
import sys
import subprocess
import os
import glob
import shutil
from pydsm._tish import _tish, _pinput
from pydsm import rootdsm_sh
import time


def test_tish():
    parameter_file = rootdsm_sh + 'AK135_SH_64.inf'
    inputs = _pinput(parameter_file)
    write_to_file = False
    start = time.time()
    u = _tish(*inputs, write_to_file)
    end = time.time()
    print('tish finished in {} s'.format(end - start))
    error_re = abs(u[2, 0, -1].real + 4.7034875e-10) / 4.7034875e-10
    error_im = abs(u[2, 0, -1].imag + 1.7664374e-11) / 1.7664374e-11
    assert error_re < 1e-7
    assert error_im < 1e-7
    print("All passed")


if __name__ == '__main__':
    test_tish()
