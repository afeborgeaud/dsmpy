import os
from dsmpy import dsm, seismicmodel
from dsmpy.event import Event
from dsmpy.station import Station
from dsmpy.utils.cmtcatalog import read_catalog
import matplotlib.pyplot as plt
import numpy as np
# load gcmt catalog
catalog = read_catalog()
event = Event.event_from_catalog(catalog,'201103171254A')

# define station FCC
stations = [
    Station(
        name='FCC', network='CN',
        latitude=58.7592, longitude=-94.0884),
    ]

plt.figure(figsize=[20, 10])

for i in range(0,9):
    os.system(f'mkdir ./data_directory/model_{i}')
    if i==0 :
        seismic_model = seismicmodel.SeismicModel.prem()
    else:
        seismic_model = seismicmodel.SeismicModel.mod_prem(i,128)    
    tlen = 3276.8 # duration of synthetics (s)
    nspc = 512 # number of points in frequency domain
    sampling_hz = 20 # sampling frequency for sythetics
    # create input parameters for pydsm
    print('seismic_model.Seismic')
    input = dsm.PyDSMInput.input_from_arrays(
        event, stations, seismic_model, tlen, nspc, sampling_hz)
    # compute synthetics in frequency domain calling DSM Fortran
    output = dsm.compute(input)
    output.to_time_domain() # perform inverse FFTAAAAAAAAAA')
    output.filter(freq=0.12) # apply a 25 seconds low-pass filter
    us = output.us # synthetics. us.shape = (3,nr,tlen)
    ts = output.ts # time points [0, tlen]
    # brackets can be used to access component and station
    u_Z_FCC = output['Z', f'FCC_CN']
    # to plot a three-component record section, use
    output.plot()
    plt.savefig(f'./data_directory/model_{i}',dpi=200)
    plt.clf()
    # to write synthetics to SAC files, use
    output.write(root_path=f'./data_directory/model_{i}',format='sac')
