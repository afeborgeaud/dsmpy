Compiled numpy.f2py (gfortran)
- numpy 1.18.1
- GNU GCC 8.3.1

How to build:

1)
```python
python -m numpy.f2py -c parameters.f90 tish.f90 others.f90 calmat.f90 trialf.f90 dclisb.f90 dclisb3.f90 -m tish
```

or

2) Run the testing script. It will first build the module.
```python
python test_tish.py
```

Usage:
```python
import tish
inputs = tish.pinput_fromfile('AK135_SH.inf')
write_to_file = False
u = tish.tish(*inputs, write_to_file) # u.shape = (3,nr,imax+1)
```
