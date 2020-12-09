from pydsm import dsm, rootdsm_psv, rootdsm_sh
from pydsm.dataset import Dataset
from pydsm.seismicmodel import SeismicModel
from pydsm.dsm import PyDSMInputFile
import os
import numpy as np
import time
from mpi4py import MPI
import glob
import sys
import pickle
import matplotlib.pyplot as plt
from obspy import read
from pydsm.spc.stf import SourceTimeFunction

def get_sac_files(root_event_folder):
    sac_files = list(glob.iglob(os.path.join(root_event_folder, '**/*Ts')))
    return sac_files

def write_outputs(outputs, root_event_folder):
    for output in outputs:
        event_id = output.event.event_id
        filename = os.path.join(root_event_folder,
                            event_id,
                            event_id + '.pkl')
        #arr = np.array([output])
        #np.save(path, arr)
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
        input_file = PyDSMInputFile(sys.argv[1])
        params = input_file.read()
        seismic_model = SeismicModel.model_from_name(params['seismic_model'])
        tlen = params['tlen']
        nspc = params['nspc']
        sampling_hz = params['sampling_hz']
        mode = params['mode']
        verbose = params['verbose']

        start_time = time.time()
        dataset = Dataset.dataset_from_sac(params['sac_files'],
                                           verbose=verbose)
        end_time = time.time()
        if verbose >= 1:
            print('Initalizing dataset finished in {} s'
                  .format(end_time - start_time))
    else:
        params = None
        dataset = None
        seismic_model = None
        tlen = None
        nspc = None
        sampling_hz = None
        mode = None
        verbose = None
    
    # log file
    log = open('log', 'w')

    # run pydsm
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen,
                                           nspc,
                                           sampling_hz,
                                           comm, mode=mode,
                                           verbose=verbose,
                                           log=log)
    end_time = time.time()
    print('rank {}: DSM finished in {} s'
          .format(rank, end_time-start_time))

    if rank == 0:
        start_time = time.time()
        for output in outputs:
            output.set_source_time_function(None)
            output.to_time_domain()
            output.write(params['output_folder'], format='sac')
        end_time = time.time()
        print('rank{}: finished FFT and writing in {} s'
              .format(rank, end_time-start_time))

        #plot(outputs[0], params['sac_files'][0])
    
    log.close()