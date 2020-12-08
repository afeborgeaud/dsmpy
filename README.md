# pydsm
Python package for computation of synthetic waveforms in spherically homogeneous transversely isotropic (VTI) media using the Direct Solution Method (DSM; [Kawai et al. 2006](https://doi.org/10.1111/j.1365-246X.2005.02829.x)).<br/><br/>
For documentation and tutorials for *pydsm*, visit the [pydsm doc page](https://afeborgeaud.github.io/pydsm/).

## INSTALLATION
1) From the Terminal, clone the pydsm git repository:
```
git clone git@github.com:afeborgeaud/pydsm.git
```

2) Set the PYTHONPATH
```export PYTHONPATH="$PYTHONPATH:<path_of_pydsm_folder>"```  
- Warning: at the moment, the pydsm folder must be named `pydsm' for the python imports to work (not, e.g., pydsm-master)

3) Install [Anaconda (python 3)](https://www.anaconda.com/products/individual)

4) Using ```conda```, install the following python packages required to run *pydsm*:
```shell
# create environment pydsm (or any other name)
conda create -n pydsm python=3.7
# install dependencies
conda install -n pydsm numpy mpi4py pandas matplotlib -y
conda install -n pydsm -c conda-forge obspy geographiclib -y
```

5) Activate the pydsm conda environment. *pydsm* should be run within this enviroment (otherwise, you might miss some dependencies)
```
conda activate pydsm
```

6) Check that pydsm has been setup succesfully:
```
python -c "import pydsm"
```
Note: Fortran sources for the DSM will be compiled the first time *pydsm* is imported. Python libraries are created from the DSM Fortran sources using numpy.f2py and the gfortran compiler. If you get compilation errors, check the following:
- gfortran >= 4.8 is required for succesful compilation, since we use the optimization flag '-Ofast'
- If you have gfortran <4.8, you should change the compiler flag from '-Ofast' to '-O3' in ```<path_of_pydsm_folder>/pydsm/__init__.py```

## GETTING STARTED
Before getting started, you should at least run ```python test_tipsv.py``` and ```python test_tish.psv``` located in in ```<path_of_pydsm_folder>/pydsm/tests```. These scripts check pydsm against pre-computed synthetics using the DSM (Fortran).

## EXAMPLES OF USAGE
1) Running pydsm using pydsm input file (run on multiple CPUs).
A template input file is in ```<path_of_pydsm_folder>/pydsm/tests/input_files/template.txt```:
```shell
sac_files ~/git/pydsm/tests/sac_files/*T
output_folder ~/git/pydsm/tests/sac_files
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
2) change the current directory to the pydsm directory
3) paste:
```shell
mpiexec -n 2 python pydsm/main.py tests/input_files/template.txt
```

2) Running pydsm from your python script.
Below is an example of python script using pydsm to compute synthetics:
```python
from pydsm import dsm, seismicmodel
from pydsm.event import Event
from pydsm.station import Station
from pydsm.utils.cmtcatalog import read_catalog
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

3) Running pydsm using a (Fortran) DSM input file.
Methods are provided to run pydsm using an input file for the (Fortran) DSM:
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
