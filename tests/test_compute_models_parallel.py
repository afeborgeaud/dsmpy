"""Test compute_models_parallel()."""

from dsmpy import dsm
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
import time
from mpi4py import MPI
import glob


def test_compute_models_parallel():
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    # Set the SAC file paths
    sac_path = "./sac_files/*T"
    sac_files = list(glob.iglob(sac_path))

    # Create the dataset
    dataset = Dataset.dataset_from_sac(sac_files, headonly=False)

    # define computation parameters
    tlen = 3276.8
    nspc = 64
    sampling_hz = 20
    mode = 2
    verbose = 2

    # Create the seismic models
    if rank == 0:
        models = [
            SeismicModel.prem(), SeismicModel.prem(), SeismicModel.ak135()]
    else:
        models = None

    # Parallel computation of synthetics
    outputs = dsm.compute_models_parallel(dataset, models,
                                          tlen, nspc, sampling_hz,
                                          comm, mode=mode, verbose=verbose)

    if rank == 0:
        assert len(outputs) == len(models)
        assert len(outputs[0]) == len(dataset.events)
        assert outputs[0][0].spcs.shape == (3, dataset.nr, nspc+1)

    return outputs, comm


if __name__ == '__main__':
    start_time = time.time()
    outputs, comm = test_compute_models_parallel()
    end_time = time.time()

    if comm.Get_rank() == 0:
        print('DSM on {} core(s) finished in {} s'
              .format(comm.Get_size(), end_time - start_time))
        print('0', outputs[0][0].spcs[2, 0, 15])
        print('1', outputs[1][0].spcs[2, 0, 15])
