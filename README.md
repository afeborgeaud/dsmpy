# pyDSM
Python package for seismic waveform computation in spherically homogeneous transversely anisotropic (VTI) media using the direct solution method (DSM).<br/>
Contains:
- Python3 wrapper for the direct solution method (DSM; Kawai et al. 2006)
- Minimal Pyhon3 implementations of key utilities of Kibrary (https://github.com/kensuke1984/Kibrary)

Compiled using numpy.f2py
- numpy 1.18.1
- GNU GCC 8.3.1

## Installation
- export PYTHONPATH=$PYTHONPATH:<path_of_pyDSM_folder> <br/>
*Warning: at the moment, the pyDSM folder must be named `pyDSM' for the imports to work (not, e.g., pyDSM-master)*

## Test
- test scripts are in pyDSM/tests

## Usage
```python
from pyDSM import dsm, rootdsm
parameter_file = rootdsm + 'AK135_SH.inf'
inputs = dsm.DSMinput(parameter_file)
outputs = dsm.compute(inputs)
outputs.to_time_domain()
u = outputs.u # u.shape = (3,nr,imax+1)
ts = outputs.ts # time points
```

Dependencies (soon):
- fftw3
