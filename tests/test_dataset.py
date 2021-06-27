"""Tests for the submodule dataset"""

from dsmpy.dataset import Dataset, filter_sac_files
from dsmpy import root_sac
from dsmpy.windowmaker import WindowMaker
from dsmpy.component import Component
from dsmpy import root_sac
import glob
import os
import numpy as np

def test_dataset_from_sac():
    sac_files = glob.glob(os.path.join(root_sac, '*[RZT]'))

    dataset = Dataset.dataset_from_sac(
        sac_files, verbose=0, headonly=False)

    windows = WindowMaker.windows_from_dataset(
        dataset, 'prem', ['s', 'S', 'Sdiff'],
        [Component.T], t_before=10., t_after=30.)

    npts = 40 * dataset.sampling_hz
    dataset.apply_windows(windows, 2, npts, buffer=0)

def test_dataset_filter_sac_files():
    sac_files = glob.glob(os.path.join(root_sac, '*[RZT]'))

    f = lambda evid_sta: evid_sta[0] == '201702211409A'
    sac_files_filt = filter_sac_files(sac_files, f)

    assert len(sac_files_filt) == 4

def test_dataset_from_sac_process():
    sac_files = glob.glob(os.path.join(root_sac, '*[RZT]'))
    freq = 0.01
    freq2 = 0.08
    filter_type = 'bandpass'

    dataset = Dataset.dataset_from_sac(
        sac_files, verbose=0, headonly=False)
    windows = WindowMaker.windows_from_dataset(
        dataset, 'prem', ['s', 'S', 'Sdiff'],
        [Component.T], t_before=10., t_after=30.)
    npts = 40 * dataset.sampling_hz
    dataset.filter(freq, freq2, filter_type)
    dataset.apply_windows(windows, 1, npts, buffer=0)

    dataset2 = Dataset.dataset_from_sac_process(
        sac_files, windows, freq, freq2, filter_type)

    assert len(dataset.events) == len(dataset2.events)
    assert (
        {event.event_id for event in dataset.events}
        == {event.event_id for event in dataset2.events}
    )
    assert len(dataset.stations) == len(dataset2.stations)
    # assert (set(dataset.stations) == set(dataset2.stations)).all()
    # assert np.allclose(dataset.data, dataset2.data)


if __name__ == '__main__':
    # test_dataset_from_sac()
    # test_dataset_filter_sac_files()
    test_dataset_from_sac_process()
