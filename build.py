import sys
import os
import subprocess
import shutil
import glob


ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dsmpy/')

def build_module(module_name='tish'):
    libroot = os.path.join(ROOT_DIR, 'flib/')
    try:
        os.mkdir(libroot)
    except FileExistsError:
        pass
    if module_name == 'tish':
        _build_tish(libroot)
    elif module_name == 'tipsv':
        _build_tipsv(libroot)
    else:
        raise RuntimeError('invalid module name')


def _build_tish(libroot):
    code = ("import sys; sys.path = {}; import numpy.f2py as f2py2e; "
            "f2py2e.main()".format(sys.path))
    sources = ['parameters.f90', 'tish.f90', 'others.f90', 'calmat.f90',
               'trialf.f90', 'dclisb.f90']
    # TODO -Ofast is a gfortran flag. Consider separate cases for other
    # compilers
    flags = ['--noopt', '--f90flags=-Ofast', '--f77flags=-Ofast']
    f2py_opts = ['-c', '-m', 'tish'] + flags + sources
    root = os.path.join(ROOT_DIR, 'src_f90/tish/')
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


def _build_tipsv(libroot):
    code = ("import sys; sys.path = {}; import numpy.f2py as f2py2e; "
            "f2py2e.main()".format(sys.path))
    sources = ['parameters.f90', 'tipsv.f90', 'others.f90', 'calmat.f90',
               'trialf.f90', 'dcsymbdl.f90', 'glu2.f90', 'rk3.f90']
    # TODO -Ofast is a gfortran flag for gcc >4.7? (not in 4.4.7).
    # Consider separate cases for older version and other compilers
    # For gcc 4.4.7 (ies clusters), use -O3 -fast-math
    flags = ['--noopt', '--f90flags=-Ofast', '--f77flags=-Ofast']
    f2py_opts = ['-c', '-m', 'tipsv'] + flags + sources
    root = os.path.join(ROOT_DIR, 'src_f90/tipsv/')
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
    import dsmpy.flib.tish
except ModuleNotFoundError:
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, '../flib')
        sys.path.append(path)
        import dsmpy.flib.tish
    except ModuleNotFoundError:
        print("Compiling tish.so from Fortran sources")
        build_module()
        try:
            import dsmpy.flib.tish
        except ModuleNotFoundError:
            raise RuntimeError("Failed to import tish.so")

try:
    import dsmpy.flib.tipsv
except ModuleNotFoundError:
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, '../flib')
        sys.path.append(path)
        import dsmpy.flib.tipsv
    except ModuleNotFoundError:
        print("Compiling tipsv.so from Fortran sources")
        build_module('tipsv')
        try:
            import dsmpy.flib.tipsv
        except ModuleNotFoundError:
            raise RuntimeError("Failed to import tish.so")