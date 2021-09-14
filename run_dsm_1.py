from dsmpy import dsm, seismicmodel
from dsmpy.event import Event
from dsmpy.station import Station
from dsmpy.utils.cmtcatalog import read_catalog
import matplotlib.pyplot as plt
# load gcmt catalog
catalog = read_catalog()
event = Event.event_from_catalog(catalog,'201103171254A')

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
