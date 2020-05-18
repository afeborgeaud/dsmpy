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
    end = time.time()
    print('tish finished in {} s'.format(end - start))

    # TODO change below for P-SV
    error_re = abs(u[2, 0, -1].real + 4.7034875e-10) / 4.7034875e-10
    error_im = abs(u[2, 0, -1].imag + 1.7664374e-11) / 1.7664374e-11
    # TODO until here

    assert error_re < 1e-7
    assert error_im < 1e-7
    print("All passed")

if __name__ == '__main__':
    test_tipsv()
