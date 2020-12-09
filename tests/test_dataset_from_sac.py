from dsmpy.dataset import Dataset
from dsmpy import root_sac
from dsmpy.windowmaker import WindowMaker
from dsmpy.component import Component
import glob
import os

if __name__ == '__main__':
    path = '/home/anselme/Dropbox/Kenji/MTZ_JAPAN/DATA/201708161251A'
    sac_files = glob.glob(os.path.join(path, '*[RZT]'))

    dataset = Dataset.dataset_from_sac(
        sac_files, verbose=2, headonly=False)

    windows = WindowMaker.windows_from_dataset(
        dataset, 'prem', ['s', 'S', 'Sdiff'],
        [Component.T], t_before=10., t_after=30.)

    npts = 40 * dataset.sampling_hz
    dataset.apply_windows(windows, 2, npts, buffer=0)
    print(dataset.data.shape)

    
