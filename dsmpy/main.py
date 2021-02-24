"""Main script for parallel synthetic computation.
The script takes a PyDSMInputFile as its single argument.

Examples:
    To run `inputfile` using 5 cores, do:
    %mpiexec -n 2 python main.py inputfile

    A template input file can be found in
    ../tests/input_files/template.txt. Run it with:
    
    %mpiexec -n 2 python main.py ../tests/input_files/template.txt

"""

from dsmpy import dsm, rootdsm_psv, rootdsm_sh
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
from dsmpy.dsm import PyDSMInputFile
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

def main(path):
    """Compute synthetics in parallel.

    Args:
        path (str): path to an input file.

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # open log file
    log = open('log', 'w', buffering=1)

    if rank == 0:
        input_file = PyDSMInputFile(path)
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
            log.write('Initalizing dataset finished in {} s\n'
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

    # run pydsm
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen,
                                           nspc,
                                           sampling_hz,
                                           mode=mode,
                                           verbose=verbose,
                                           log=log)
    end_time = time.time()
    log.write('rank {}: DSM finished in {} s\n'
          .format(rank, end_time-start_time))

    if rank == 0:
        start_time = time.time()
        for output in outputs:
            output.set_source_time_function(None)
            output.to_time_domain()
            output.write(params['output_folder'], format='sac')
        end_time = time.time()
        log.write('finished FFT and writing in {} s\n'
              .format(rank, end_time-start_time))
    
    log.close()
    
    return "Done!"

if __name__ == '__main__':
    path = sys.argv[1]
    main(path)