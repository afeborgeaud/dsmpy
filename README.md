# dsmpy
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![PyPI version fury.io](https://d25lcipzij17d.cloudfront.net/badge.svg?id=py&type=6&v=1.0a4&x2=0)](https://test.pypi.org/project/dsmpy/)

Python package for computation of synthetic waveforms in spherically homogeneous transversely isotropic (VTI) media using the Direct Solution Method (DSM; [Kawai et al. 2006](https://doi.org/10.1111/j.1365-246X.2005.02829.x)).<br/>
The original DSM softwares can be found on the [Github page of UT Global Seismology](https://github.com/UT-GlobalSeismology).<br/><br/>
The documentation for dsmpy with several usage examples can be found [here](https://afeborgeaud.github.io/dsmpy/).

# INSTALLATION
## Preferred method: dependencies using conda and building from source
1) Clone the dsmpy repository
```bash
git clone https://github.com/seismobassoon/dsmpy
```
2) Install conda-build in ```base```
```bash
conda install conda-build -y
```
3) Create the ```dsm``` conda environment using the ```environment.yml``` YAML file:
```bash
conda env create -f /path/to/dsmpy/environment.yml
```
3bis) You can add the yml file to your environment by:
```bash
conda env update --prefix ./env --file environment.yml  --prune
```
Anyways, be careful with the compatibility of fortran and c(lang) compilers: it's secured to install them via:
```
conda install -c conda-forge c-compiler
```
4) Build dsmpy. ```/path/to/dsmpy/``` is the path to the local dsmpy git repository:
```bash
conda develop -n dsm -b /path/to/dsmpy/
```
4b) OK if your conda cannot do this by reasons XYZ, try this:
```python setup.py install```
5) Run tests. From /path/to/dsmpy/
```bash
conda activate dsm

6) Modified (might not be necessay for all, but for some systems): In case of the module absence error (tish, tipsv) during pytest
python build.py

pytest
```

## Using conda to install dsmpy conda package (currently not the latest version)
We recommend using ```conda``` to install ```dsmpy```, as it takes care of the dependencies required to compile the Fortran sources in dsmpy. At the moment, dsmpy has been compiled for ```linux-64``` and ```osx-64``` platforms. <br>
To install ```dsmpy``` using ```conda```:
1) From a terminal, run
```bash
conda create -n dsmpy
conda install -n dsmpy -c afeborgeaud -c conda-forge dsmpy
```

## Using pip (TestPyPI)
1) (Optional) You may want to install dsmpy in a virtual environment. If so, do
```
python3 -m venv venv
source ./venv/bin/activate
```
2) Install [*build*](https://pypi.org/project/build/), a PEP517 package builder
```
pip install build
```
3) In a shell, type
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple dsmpy
```
This will download the dsmpy package from the TestPyPI respository and the required dependencies from the PyPI repository.

4) Check that dsmpy has been installed succesfully:
```
python -c "import dsmpy"
```
**Note:** Fortran sources for the DSM will be compiled during the installation (using numpy.f2py and the GNU Fortran compiler). If you get compilation errors, check the following:
- gfortran>=4.8 is required for succesful compilation, because of the optimization flag '-Ofast'
- If you have gfortran<4.8, you should change the compiler flag from '-Ofast' to '-O3' in ```<path_of_dsmpy_folder>/pydsm/__init__.py```

## Build from source using pip
1) Clone the dsmpy repository
```
git clone https://github.com/afeborgeaud/dsmpy
```
2) (Optional) You may want to install dsmpy in a virtual environment. If so, do
```
python3 -m venv venv
source ./venv/bin/activate
```
3) Install requirements
```
python3 -m pip install -r requirements.txt
```
5) Install [*build*](https://pypi.org/project/build/), a PEP517 package builder
```
pip install build
```
4) To build the dsmpy package, from the root directory ```dsmpy``` run
```
python -m build .
```
5) This creates ```.whl``` and ```.gz.tar``` dist files in the ```dist``` directory. Now pydsm can be installed with
```
pip install dist/*.whl
```
or
```
pip install dist/*.tar.gz
```

# EXAMPLES
1) Running dsmpy using an input file (run on multiple CPUs).
A template input file is in ```<path_of_pydsm_folder>/dsmpy/tests/input_files/template.txt```:
```shell
sac_files ~/git/dsmpy/tests/sac_files/*T
output_folder ~/git/dsmpy/tests/sac_files
# duration of synthetics (in seconds)
tlen 3276.8
# number of points of frequency-domain synthetics
# minimum period Tmin = tlen / nspc (s)
nspc 256 
# sampling frequency for time-domain synthetics
sampling_hz 20
# prem, ak135
seismic_model prem 
# 0: P-SV+SH, 1: P-SV, 2: SH (default: 0)
mode 0
# 0: quiet, 1: talkative, 2: debug (default: 0)
verbose 0
```

To run this input file on 2 CPUs:
1) open a Terminal 
2) change the current directory to the dsmpy directory
3) paste:
```shell
mpiexec -n 2 python pydsm/main.py tests/input_files/template.txt
```

2) Running dsmpy from a python script.
Below is an example of python script using dsmpy to compute synthetics:
```python
from dsmpy import dsm, seismicmodel
from dsmpy.event import Event
from dsmpy.station import Station
from dsmpy.utils.cmtcatalog import read_catalog
# load gcmt catalog
catalog = read_catalog()
# get event from catalog
event = Event.event_from_catalog(
    catalog, '200707211534A')
# define station FCC
stations = [
    Station(
        name='FCC', network='CN',
        latitude=58.7592, longitude=-94.0884), 
    ]
# load (anisotropic) PREM model
seismic_model = seismicmodel.SeismicModel.prem()
tlen = 3276.8 # duration of synthetics (s)
nspc = 256 # number of points in frequency domain
sampling_hz = 20 # sampling frequency for sythetics
# create input parameters for pydsm
input = dsm.PyDSMInput.input_from_arrays(
    event, stations, seismic_model, tlen, nspc, sampling_hz)
# compute synthetics in frequency domain calling DSM Fortran
output = dsm.compute(input)
output.to_time_domain() # perform inverse FFT
output.filter(freq=0.04) # apply a 25 seconds low-pass filter
us = output.us # synthetics. us.shape = (3,nr,tlen)
ts = output.ts # time points [0, tlen]
# brackets can be used to access component and station
u_Z_FCC = output['Z', 'FCC_CN']
# to plot a three-component record section, use
output.plot()
plt.show()
# to write synthetics to SAC files, use
output.write(root_path='.', format='sac')
```

3) Running dsmpy using a (Fortran) DSM input file.
Methods are provided to run dsmpy using an input file for the (Fortran) DSM:
```python
from pydsm import dsm, rootdsm_sh
parameter_file = rootdsm_sh + 'AK135_SH.inf'
inputs = dsm.PyDSMInput.input_from_file(parameter_file, mode=2)
outputs = dsm.compute(inputs, mode=2)
outputs.to_time_domain()
us = outputs.us    # us.shape = (3,nr,tlen)
ts = outputs.ts    # len(ts) = tlen
stations = outputs.stations        # len(stations) = nr
components = outputs.components    # len(components) = 3
```
