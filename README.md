Compiled with:

python -m numpy.f2py -c parameters.f90 tish.f90 others.f90 calmat.f90 trialf.f90 dclisb.f90 dclisb3.f90 -m tish

>>> import tish
>>> inputs = tish.pinput_fromfile('AK135_SH.inf')
>>> write_to_file = False
>>> u = tish.tish(*inputs, write_to_file) # u.shape = (3,nr,imax+1)
