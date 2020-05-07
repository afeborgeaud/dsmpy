import numpy as np
from numpy import f2py
import sys
import subprocess
import os
import glob
import shutil
from pydsm._tish import _tish, _pinput

def test_tish():
    inputs = _pinput(
        '../pyDSM/src_f90/tish/example/dsm_accuracy_check/AK135_SH_64.inf')
    write_to_file = False
    u = _tish(*inputs, write_to_file)
    error_re = abs(u[2, 0, -1].real + 4.7034875e-10) / 4.7034875e-10
    error_im = abs(u[2, 0, -1].imag + 1.7664374e-11) / 1.7664374e-11
    assert error_re < 1e-7
    assert error_im < 1e-7
    print("All passed")

if __name__ == '__main__':
    test_tish()