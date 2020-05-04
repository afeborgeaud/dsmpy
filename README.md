Compiled numpy.f2py (gfortran)
- numpy 1.18.1
- GNU GCC 8.3.1

**Setup:**
- export PYTHONPATH=$PYTHONPATH:<path_of_pyDSM_folder>

**Tests:**
- test scripts are in pyDSM/tests

**Usage:**
```python
from pyDSM import dsm
inputs = dsm.DSMinput(parameter_file)
outputs = dsm.compute(inputs)
outputs.to_time_domain()
u = outputs.u # u.shape = (3,nr,imax+1)
ts = outputs.ts # time points
```

Dependencies (soon):
- fftw3
