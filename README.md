# dsmpy
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![PyPI version fury.io](https://d25lcipzij17d.cloudfront.net/badge.svg?id=py&type=6&v=1.0a4&x2=0)](https://test.pypi.org/project/dsmpy/)

Python package for computation of synthetic waveforms in spherically homogeneous transversely isotropic (VTI) media using the Direct Solution Method (DSM; [Kawai et al. 2006](https://doi.org/10.1111/j.1365-246X.2005.02829.x)).<br/>
The original DSM softwares can be found on the [Github page of UT Global Seismology](https://github.com/UT-GlobalSeismology).<br/><br/>
The documentation for dsmpy with several usage examples can be found [here](https://afeborgeaud.github.io/dsmpy/).

# For Linux users
At installation, dsmpy needs to compile Fortran and C librairies, and needs openmpi for parallel computing.
Currently, dsmpy has been tested with gcc.

# For Windows users
At installation, dsmpy needs to compile Fortran and C librairies, and needs openmpi for parallel computing.
The installation on Windows machines has been tested on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install])

A quick summary to set up the environment tested:

1) Open a PowerShell (in admin mode)

2) In the PowerShell, type
```
wsl --install
```
This should install a ubuntu distribution by default, and launch a Ubuntu bash terminal 

3) From the Ubuntu terminal, install python, gcc and openmpi
```
sudo apt-get update && apt-get install -y python3 python3-pip
sudo apt install python-is-python3
sudo apt-get install gcc
sudo apt-get install -y openmpi-bin libopenmpi-dev
```

4) Create a directory for your python project (here we assume you have a ```git``` folder in your home directory), and open it in Visual Studio Code
```
cd ~/git
mkdir my_project
# open my_project in VS Code
code my_project
```
Your are ready to proceed to the installation of dsmpy


# INSTALLATION
## Prefererd method: install using pip from github
1) (Optional) You may want to install dsmpy in the virtual environment [env_name]. If so, from your project in a terminal:
```
python3 -m venv .venvs/[env_name]
source .venvs/bin/[env_name]/activate
```
2)
```
pip install dsmpy@git+https://github.com/afeborgeaud/dsmpy@v1.0a5
```

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
mpiexec -n 2 python dsmpy/main.py tests/input_files/template.txt
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
plt.savefig('output.png')
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
