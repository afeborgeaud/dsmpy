from pydsm import dsm, rootdsm_psv
from pydsm.dataset import Dataset
from pydsm.seismicmodel import SeismicModel
import os
import numpy as np
import time
from mpi4py import MPI
import glob
import sys

def run(comm, dataset, seismic_):
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        parameter_files = [
        rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf']
        dataset = Dataset.dataset_from_files(parameter_files)
        seismic_model = SeismicModel.prem()
        tlen = 3276.8
        nspc = 64
        sampling_hz = 20
    else:
        dataset = None
        seismic_model = None
        tlen = None
        nspc = None
        sampling_hz = None
    
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen, nspc, sampling_hz,
                                           comm)
    end_time = time.time()

    return outputs

def get_sac_files(root_event_folder):
    sac_files = list(glob.iglob(os.path.join(root_event_folder, '**/*T')))
    return sac_files

def write_outputs(outputs, root_event_folder):
    for output in outputs:
        event_id = outputs.event.event_id
        path = os.path.join(root_event_folder,
                            event_id,
                            event_id)
        arr = np.array(output)
        np.write(path, arr)

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        #root_event_folder = sys.argv[1]
        root_event_folder = '/work/anselme/TEST_PYDSM/test'
        sac_files = get_sac_files(root_event_folder)

    if rank == 0:
        dataset = Dataset.dataset_from_sac(sac_files)
        seismic_model = SeismicModel.prem()
        tlen = 3276.8
        nspc = 64
        sampling_hz = 20
    else:
        dataset = None
        seismic_model = None
        tlen = None
        nspc = None
        sampling_hz = None
    
    # run DSM
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen, nspc, sampling_hz,
                                           comm)
    end_time = time.time()
    print('rank {}: DSM finished in {} s'
          .format(rank, end_time-start_time))
    
    if rank == 0:
        write_outputs(outputs, root_event_folder)
    