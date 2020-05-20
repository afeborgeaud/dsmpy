# pydsm
Python package for parallel seismic waveform computation in spherically homogeneous transversely anisotropic (VTI) media using the direct solution method (DSM).<br/>
Contains:
- Python3 wrapper for the direct solution method (DSM; Kawai et al. 2006)
- Minimal Python3 implementations of key utilities of Kibrary (https://github.com/kensuke1984/Kibrary)

Compiled using numpy.f2py
- numpy 1.18.1
- GNU GCC 8.3.1

## Installation
- export PYTHONPATH=$PYTHONPATH:<path_of_pydsm_folder> <br/>
*Warning: at the moment, the pydsm folder must be named `pydsm' for the imports to work (not, e.g., pydsm-master)*

## Test
- test scripts are in pydsm/tests

## Usage
#### Using pydsm Python classes
```python
from pydsm import dsm, seismicmodel
from pydsm.utils.cmtcatalog import read_catalog
# load gcmt catalog
catalog = read_catalog()
# get event from catalog
event = dsm.Event.event_from_catalog(
    catalog, '200707211534A')
# define station FCC
stations = [
    dsm.Station(
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
```

#### Using Fortran DSM input file
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

## Dependencies:
- numpy
- mpi4py
- pandas
- matplotlib
- obspy
- fftw3 (soon)
