import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name="dsmpy",
    version='2020.12',
    author='Anselme Borgeaud, Kensuke Konishi',
    author_email='aborgeaud@gmail.com',
    description='Python wrapper for DSM',
    long_description = long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'obspy',
        'numpy',
        'mpi4py',
        'pandas',
        'matplotlib',
        'geographiclib',
    ],
    python_requires='>=3.7'
)


from numpy.distutils.core import setup, Extension

lib_tish = Extension(
    name='pydsm.flib.tish',
    sources=[
        'pydsm/src_f90/tish/parameters.f90',
        'pydsm/src_f90/tish/tish.f90',
        'pydsm/src_f90/tish/others.f90',
        'pydsm/src_f90/tish/calmat.f90',
        'pydsm/src_f90/tish/trialf.f90',
        'pydsm/src_f90/tish/dclisb.f90'
    ],
    extra_f90_compile_args=['-Ofast'],
    extra_f77_compile_args=['-Ofast'],
)

lib_tipsv = Extension(
    name='pydsm.flib.tipsv',
    sources=[
        'pydsm/src_f90/tipsv/parameters.f90',
        'pydsm/src_f90/tipsv/tipsv.f90',
        'pydsm/src_f90/tipsv/others.f90',
        'pydsm/src_f90/tipsv/calmat.f90',
        'pydsm/src_f90/tipsv/trialf.f90',
        'pydsm/src_f90/tipsv/dcsymbdl.f90',
        'pydsm/src_f90/tipsv/glu2.f90',
        'pydsm/src_f90/tipsv/rk3.f90'
    ],
    extra_f90_compile_args=['-Ofast'],
    extra_f77_compile_args=['-Ofast'],
)

setup(
    ext_modules=[lib_tish, lib_tipsv],
)