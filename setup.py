import setuptools
from numpy.distutils.core import setup, Extension
from distutils.command.sdist import sdist

lib_tish = Extension(
        name='dsmpy.flib.tish',
        sources=[
            'dsmpy/src_f90/tish/parameters.f90',
            'dsmpy/src_f90/tish/tish.f90',
            'dsmpy/src_f90/tish/others.f90',
            'dsmpy/src_f90/tish/calmat.f90',
            'dsmpy/src_f90/tish/trialf.f90',
            'dsmpy/src_f90/tish/dclisb.f90'
        ],
        extra_f90_compile_args=['-Ofast'],
        extra_f77_compile_args=['-Ofast'],
    )

lib_tipsv = Extension(
    name='dsmpy.flib.tipsv',
    sources=[
        'dsmpy/src_f90/tipsv/parameters.f90',
        'dsmpy/src_f90/tipsv/tipsv.f90',
        'dsmpy/src_f90/tipsv/others.f90',
        'dsmpy/src_f90/tipsv/calmat.f90',
        'dsmpy/src_f90/tipsv/trialf.f90',
        'dsmpy/src_f90/tipsv/dcsymbdl.f90',
        'dsmpy/src_f90/tipsv/glu2.f90',
        'dsmpy/src_f90/tipsv/rk3.f90'
    ],
    extra_f90_compile_args=['-Ofast'],
    extra_f77_compile_args=['-Ofast'],
)

if __name__ == '__main__':
    with open('README.md', 'r') as fh:
        long_description = fh.read()

    setup(
        name="dsmpy",
        version='1.0a4',
        author='Anselme Borgeaud, Kensuke Konishi, Nobuaki Fuji',
        author_email='nobuaki@ipgp.fr',
        license='MIT',
        description='Python wrapper for DSM',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/afeborgeaud/dsmpy',
        packages=setuptools.find_packages(),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Developers",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        install_requires=[
            'obspy',
            'mpi4py',
            'pandas',
            'matplotlib',
            'geographiclib',
            "numpy",
            "pytest",
        ],
        ext_modules=[lib_tish, lib_tipsv],
        python_requires='>=3.9',
        package_data={
            'dsmpy' : ['resources/scardec.pkl'],
        },
        cmdclass={'sdist': sdist},
    )
