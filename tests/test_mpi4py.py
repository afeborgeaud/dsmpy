from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_cores = comm.Get_size()

if rank == 0:
    arr = np.arange(24, dtype=np.float64)
    sendcounts = (6, 6, 12)
    displacements = (0, 6, 12)
    arr_local = np.empty(6)
else:
    arr = None
    sendcounts = None
    displacements = None
    if rank == 1:
        arr_local = np.empty(6)
    else:
        arr_local = np.empty(12)

comm.Scatterv([arr, sendcounts, displacements, MPI.DOUBLE],
              arr_local, root=0)

if rank == 2:
    arr_local = arr_local.reshape((2, 2, 3))
else:
    arr_local = arr_local.reshape((2, 1, 3))

comm.Barrier()

if rank == 0:
    arr_gathered = np.empty((2, 4, 3), dtype=np.float64)
else:
    arr_gathered = None

comm.Gatherv(arr_local,
            [arr_gathered, sendcounts, displacements, MPI.DOUBLE],
            root=0)

if rank == 0:
    print(arr_gathered)