from pydsm import dsm, rootdsm_psv
import os
import numpy as np
import time
from mpi4py import MPI

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        parameter_file = os.path.join(rootdsm_psv, 'test3.inf')
        print(parameter_file)
        inputs = dsm.PyDSMInput.input_from_file(parameter_file, mode=1)
    else:
        inputs = None

    start_time = time.time()
    spcs = dsm.compute_parallel(inputs, comm)
    end_time = time.time()

    if rank == 0:
        print('DSM on {} cores finished in {} s'
              .format(n_cores, end_time - start_time))
    
    if rank == 0 and n_cores == 1:
        np.save('spcs_n1', spcs)
    elif rank == 0:
        spcs_n1 = np.load('spcs_n1.npy')
        assert np.allclose(spcs_n1, spcs, rtol=1e-17)
        print("Same!")
