from dsmpy import dsm, rootdsm_psv
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
from dsmpy.station import Station
from dsmpy.event import Event
from dsmpy.utils.cmtcatalog import read_catalog
from dsmpy._tish import parameters
import os
import numpy as np
import time
from mpi4py import MPI

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    tlen = 3276.8
    nspc = 32
    sampling_hz = 20
    maxnzone = parameters['maxnzone']

    if rank == 0:
        models = [SeismicModel.prem(), SeismicModel.ak135()]
    else:
        models = None

    catalog = read_catalog()
    event = Event.event_from_catalog(
        catalog, '200707211534A')
    stations = [
        Station(
            name='FCC', network='CN',
            latitude=58.7592, longitude=-94.0884), 
        ]
    events = [event]
    stations = [stations]
    dataset = Dataset.dataset_from_arrays(
        events, stations, sampling_hz)
    
    start_time = time.time()
    outputs = dsm.compute_models_parallel(dataset, models,
                                           tlen, nspc, sampling_hz,
                                           comm, maxnzone, mode=0, verbose=2)
    end_time = time.time()

    if rank == 0:
        print('DSM on {} cores finished in {} s'
              .format(n_cores, end_time - start_time))
        
    if rank == 0 :
        print('0', outputs[0][0].spcs[2, 15, 0])
        print('1', outputs[1][0].spcs[2, 15, 0])
