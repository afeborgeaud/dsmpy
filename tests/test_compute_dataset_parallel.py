from dsmpy import dsm, rootdsm_psv
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
from dsmpy.component import Component
from dsmpy import root_sac
import os
import numpy as np
import time
from mpi4py import MPI
import matplotlib.pyplot as plt
import glob

def test_compute_models_parallel():
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    # Set the SAC file paths
    sac_files = list(glob.iglob(os.path.join(root_sac, '*2017*T')))

    # Create the dataset
    dataset = Dataset.dataset_from_sac(sac_files, headonly=False)

    # define computation parameters
    tlen = 1638.4
    nspc = 64
    sampling_hz = 20
    mode = 2

    models = [SeismicModel.prem()]

    # Parallel computation of synthetics
    outputs_1 = dsm.compute_models_parallel(
        dataset, models, tlen, nspc, sampling_hz, mode)

    outputs_2 = dsm.compute_dataset_parallel(
        dataset, models[0], tlen, nspc, sampling_hz, mode)

    if rank == 0:
        assert np.allclose(outputs_1[0][0].spcs, outputs_2[0].spcs)

if __name__ == '__main__':
    test_compute_models_parallel()