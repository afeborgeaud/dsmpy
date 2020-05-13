from pydsm import dsm, rootdsm_psv
import os
import numpy as np
import time

if __name__ == '__main__':
    parameter_file = os.path.join(rootdsm_psv, 'test2.inf')
    inputs = dsm.PyDSMInput.input_from_file(parameter_file)

    start_time = time.time()
    spcs, comm = dsm.compute_parallel(inputs)
    end_time = time.time()
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == 0:
        print('DSM on {} cores finished in {} s'
            .format(size, end_time - start_time))
    #if rank == 0:
        #np.save('spcs_n1', spcs)
        #spcs_n1 = np.load('spcs_n1.npy')
        #assert np.allclose(spcs_n1, spcs, rtol=1e-17)
        #print("Same!")
