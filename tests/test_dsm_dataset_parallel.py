from pydsm import dsm, rootdsm_psv
from pydsm.dataset import Dataset
import os
import numpy as np
import time
from mpi4py import MPI

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        parameter_files = [
        rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf']
        dataset = Dataset.dataset_from_files(parameter_files)
    else:
        dataset = None
    
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, comm)
    end_time = time.time()

    if rank == 0:
        print('DSM on {} cores finished in {} s'
              .format(n_cores, end_time - start_time))
        
    if rank == 0 :
        for i, output in enumerate(outputs):
            filename = '{}_{}_n{}'.format(output.event, i, n_cores)
            np.save(filename, output.spcs)
    #elif rank == 0:
            if i == 0:
                spcs_n1 = np.load('spcs_n1.npy')
                assert np.allclose(spcs_n1, output.spcs, rtol=1e-17)
                print("Same!")
