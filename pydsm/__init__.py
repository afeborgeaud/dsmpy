import sys
import os
import subprocess
import shutil
import glob
from pydsm.definitions import ROOT_DIR


def build_module(module_name='tish'):
    if module_name == 'tish':
        _build_tish()
    elif module_name == 'tipsv':
        _build_tipsv()
    else:
        raise RuntimeError('invalid module name')


def _build_tish():
    code = ("import sys; sys.path = {}; import numpy.f2py as f2py2e; "
            "f2py2e.main()".format(sys.path))
    sources = ['parameters.f90', 'tish.f90',
               'others.f90', 'calmat.f90', 'trialf.f90',
               'dclisb.f90', 'dclisb3.f90']
    f2py_opts = ['-c', '-m', 'tish'] + sources
    root = os.path.join(ROOT_DIR, 'src_f90/tish/')
    libroot = os.path.join(ROOT_DIR, 'lib/')
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


def _build_tipsv():
    code = ("import sys; sys.path = {}; import numpy.f2py as f2py2e; "
            "f2py2e.main()".format(sys.path))
    sources = ['parameters.f90', 'tipsv.f90',
               'others.f90', 'calmat.f90', 'trialf.f90',
               'dcsymbdl.f90', 'glu2.f90', 'rk3.f90']
    f2py_opts = ['-c', '-m', 'tipsv'] + sources
    root = os.path.join(ROOT_DIR, 'src_f90/tipsv/')
    libroot = os.path.join(ROOT_DIR, 'lib/')
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
        for lib in glob.glob('tipsv.cpython-3*'):
            shutil.move(lib, os.path.join(libroot, lib))
    finally:
        os.chdir(cwd)


try:
    import pydsm.lib.tish
except ModuleNotFoundError:
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, '../lib')
        sys.path.append(path)
        import pydsm.lib.tish
    except ModuleNotFoundError:
        print("Compiling tish.so from Fortran sources")
        build_module()
        try:
            import pydsm.lib.tish
        except ModuleNotFoundError:
            raise RuntimeError("Failed to import tish.so")

try:
    import pydsm.lib.tipsv
except ModuleNotFoundError:
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, '../lib')
        sys.path.append(path)
        import pydsm.lib.tipsv
    except ModuleNotFoundError:
        print("Compiling tipsv.so from Fortran sources")
        build_module('tipsv')
        try:
            import pydsm.lib.tipsv
        except ModuleNotFoundError:
            raise RuntimeError("Failed to import tish.so")

#
rootdsm = os.path.join(ROOT_DIR,
                       'src_f90/tish/example/dsm_accuracy_check/')
