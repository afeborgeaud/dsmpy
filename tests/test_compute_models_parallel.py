"""Test compute_models_parallel()."""

from dsmpy import dsm
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
import time
from mpi4py import MPI
import glob
import numpy as np

MAX_SPC = 1e6

def test_compute_models_parallel():
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    # Set the SAC file paths
    # sac_path = "./sac_files/*T"
    sac_path = "/work/anselme/japan/DATA/200502*/*T"
    sac_files = list(glob.iglob(sac_path))

    # Create the dataset
    dataset = Dataset.dataset_from_sac(sac_files, headonly=False)

    # define computation parameters
    tlen = 3276.8
    nspc = 64
    sampling_hz = 20
    mode = 2
    verbose = 1

    # Create the seismic models
    if rank == 0:
        models = [
            SeismicModel.prem(), SeismicModel.prem(), SeismicModel.ak135()]
    else:
        models = None

    # Parallel computation of synthetics
    outputs = dsm.compute_models_parallel(dataset, models,
                                          tlen, nspc, sampling_hz,
                                          mode=mode, verbose=verbose)

    if rank == 0:
        assert len(outputs) == len(models)
        for imod in range(len(models)):
            assert len(outputs[imod]) == len(dataset.events)
            for iev in range(len(dataset.events)):
                assert outputs[imod][iev].spcs.shape == (
                    3, dataset.nrs[iev], nspc + 1)
                assert not np.isnan(outputs[imod][iev].spcs[2]).any()
                assert (
                        np.abs(outputs[imod][iev].spcs[2].max(axis=1))
                        < MAX_SPC).all()
                assert (
                    np.abs(outputs[imod][iev].spcs[2]).max(axis=1)
                    > 0).all()

    return outputs


if __name__ == '__main__':
    start_time = time.time()
    outputs, comm = test_compute_models_parallel()
    end_time = time.time()

    if MPI.COMM_WORLD.Get_rank() == 0:
        print('DSM on {} core(s) finished in {} s'
              .format(comm.Get_size(), end_time - start_time))
        print('0', outputs[0][0].spcs[2, 0, 15])
        print('1', outputs[1][0].spcs[2, 0, 15])
