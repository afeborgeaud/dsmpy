import numpy as np
from numpy import f2py
import sys
import subprocess
import os
import glob
import shutil

def build_module(module_name='tish'):
    code = ("import sys; sys.path = {}; import numpy.f2py as f2py2e; " 
            "f2py2e.main()".format(sys.path))
    sources = ['parameters.f90', 'tish.f90',
        'others.f90', 'calmat.f90', 'trialf.f90',
        'dclisb.f90', 'dclisb3.f90']
    f2py_opts = ['-c', '-m', module_name] + sources
    root = os.path.abspath('../src_f90/tish/')
    libroot = os.path.abspath('../lib/')
    cwd = os.getcwd()

    try:
        os.chdir(root)
        if 'parameters.mod' in os.listdir('.'):
            os.remove('parameters.mod')
        cmd = [sys.executable, '-c', code] + f2py_opts
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        out, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("Failed running f2py: {}\n{}"
            .format(cmd[4:], out.decode('utf-8')))
        for lib in glob.glob('tish.cpython-3*'):
            shutil.move(lib, os.path.join(libroot, lib))
    finally:
        os.chdir(cwd)

try:
    import tish
except ModuleNotFoundError:
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, '../lib')
        sys.path.append(path)
        import tish
    except ModuleNotFoundError:
        print("Compiling tish.so from Fortran sources")
        build_module()
        try:
            import tish
        except ModuleNotFoundError:
            raise RuntimeError("Failed to import tish.so")

def test_tish():
    inputs = tish.pinput_fromfile(
        '../src_f90/tish/example/dsm_accuracy_check/AK135_SH_64.inf')
    write_to_file = False
    u = tish.tish(*inputs, write_to_file)
    error_re = abs(u[2, 0, -1].real + 4.7034875e-10) / 4.7034875e-10
    error_im = abs(u[2, 0, -1].imag + 1.7664374e-11) / 1.7664374e-11
    assert error_re < 1e-7
    assert error_im < 1e-7
    print("All passed")

if __name__ == '__main__':
    test_tish()