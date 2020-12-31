"""Tests for the submodule dataset"""

from dsmpy.dataset import Dataset, filter_sac_files
from dsmpy import root_sac
from dsmpy.windowmaker import WindowMaker
from dsmpy.component import Component
import glob
import os

sac_file_root = './sac_files'

def test_dataset_from_sac():
    sac_files = glob.glob(os.path.join(sac_file_root, '*[RZT]'))

    dataset = Dataset.dataset_from_sac(
        sac_files, verbose=0, headonly=False)

    windows = WindowMaker.windows_from_dataset(
        dataset, 'prem', ['s', 'S', 'Sdiff'],
        [Component.T], t_before=10., t_after=30.)

    npts = 40 * dataset.sampling_hz
    dataset.apply_windows(windows, 2, npts, buffer=0)


def test_dataset_filter_sac_files():
    sac_files = glob.glob(os.path.join(sac_file_root, '*[RZT]'))

    f = lambda evid_sta: evid_sta[0] == '201702211409A'
    sac_files_filt = filter_sac_files(sac_files, f)

    assert len(sac_files_filt) == 4


if __name__ == '__main__':
    test_dataset_from_sac()
    test_dataset_filter_sac_files()
