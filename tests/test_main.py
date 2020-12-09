from dsmpy import dsm, rootdsm_psv, rootdsm_sh
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
from dsmpy.spc.stf import SourceTimeFunction
import os
import numpy as np
import time
from mpi4py import MPI
import glob
import sys
import pickle
import matplotlib.pyplot as plt
from obspy import read

def get_sac_files(root_event_folder):
    sac_files = list(glob.iglob(os.path.join(root_event_folder, '**/*Ts')))
    return sac_files

def write_outputs(outputs, root_event_folder):
    for output in outputs:
        event_id = output.event.event_id
        filename = os.path.join(root_event_folder,
                            event_id,
                            event_id + '.pkl')
        with open(filename, 'wb') as f:
            pickle.dump(output, f)

def plot(output, sac_file):
    trace = read(sac_file)[0]
    residual = output.us[2, 0] - trace
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)
    ax0.plot(output.ts, output.us[2, 0], color='blue', label='pydsm')
    ax0.plot(output.ts, trace.data, color='red', label='DSM')
    ax1.plot(output.ts, residual, color='blue', label='residual')
    ax1.set(xlabel='Time (s)',
           ylabel='Velocity ()')
    ax0.set(ylabel='Velocity ()')
    ax0.legend()
    ax1.legend()
    plt.show()

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        #root_event_folder = sys.argv[1]
        root_event_folder = os.path.join(rootdsm_sh, 'sac')
        sac_files = get_sac_files(root_event_folder)
        print(root_event_folder)

    if rank == 0:
        dataset = Dataset.dataset_from_sac(sac_files)
        seismic_model = SeismicModel.ak135()
        tlen = 3276.8
        nspc = 512
        sampling_hz = 20
    else:
        dataset = None
        seismic_model = None
        tlen = None
        nspc = None
        sampling_hz = None
    
    # run pydsm
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen, nspc, sampling_hz,
                                           comm)
    end_time = time.time()
    print('rank {}: DSM finished in {} s'
          .format(rank, end_time-start_time))

    for output in outputs:
        # output.set_source_time_function(
        #     SourceTimeFunction.triangle(10., tlen, nspc))
        output.to_time_domain()
    
    # debug for station order
    for station, sac in zip(
            np.concatenate([output.stations for output in outputs]),
            sac_files):
        st = read(sac)
        assert str(station) == (st[0].stats.station + '_' 
                               + st[0].stats.network)
    
    if rank == 0:
        write_outputs(outputs, root_event_folder)

    if rank == 0:
        plot(outputs[0], sac_files[0])
    