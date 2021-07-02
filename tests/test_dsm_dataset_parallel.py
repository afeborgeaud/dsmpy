from dsmpy import dsm, rootdsm_psv
from dsmpy.dataset import Dataset
from dsmpy.seismicmodel import SeismicModel
from dsmpy.component import Component
import os
import numpy as np
import time
from mpi4py import MPI
import matplotlib.pyplot as plt
import glob

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    n_cores = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        parameter_files = [
        # rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf']
        dataset = Dataset.dataset_from_files(parameter_files)
        seismic_model = SeismicModel.prem()
        tlen = 3276.8
        nspc = 64
        sampling_hz = 20

        sac_files = glob.glob(
            rootdsm_psv + '200901010139A/sac/200901010139A/*Zs')
        dataset2 = Dataset.dataset_from_sac(sac_files)
        dataset2.mts[0] = dataset.mts[0]
        dataset2.events[0].source_time_function = None
    else:
        dataset = None
        seismic_model = None
        tlen = None
        nspc = None
        sampling_hz = None
        
        dataset2 = None

    if rank == 0:
        assert np.allclose(dataset.lats, dataset2.lats)
        assert np.allclose(np.sort(dataset.lons), np.sort(dataset2.lons))
        assert np.allclose(np.sort(dataset.phis), np.sort(dataset2.phis))
        assert np.allclose(np.sort(dataset.thetas), np.sort(dataset2.thetas))
        print('Datasets identical')

        print(dataset2.stations[10].latitude, dataset2.stations[10].longitude)
        print(dataset2.lats[10], dataset2.lons[10])
    
    start_time = time.time()
    outputs = dsm.compute_dataset_parallel(dataset, seismic_model,
                                           tlen, nspc, sampling_hz,
                                           mode=2)

    outputs2 = dsm.compute_dataset_parallel(dataset2, seismic_model,
                                           tlen, nspc, sampling_hz,
                                           mode=2)
    end_time = time.time()

    if rank == 0:
        print('DSM on {} cores finished in {} s'
              .format(n_cores, end_time - start_time))

    if rank == 0:
        fig, ax = outputs2[0].plot_component(Component.T, color='black')
        # plt.show()
        _, ax = outputs[0].plot_component(Component.T, ax=ax, color='red')
        plt.show()

        for i, output in enumerate(outputs):
            filename = '{}_{}_n{}_tmp'.format(output.event, i, n_cores)
            np.save(filename, output.spcs)
    #elif rank == 0:
            if i == 1:
                spcs_n2 = np.load('20090101_1_n2_tmp.npy')
                assert np.allclose(spcs_n2, output.spcs, atol=1e-17)
                print("Same!")

            # output.plot_component(Component.T)
            # plt.show()