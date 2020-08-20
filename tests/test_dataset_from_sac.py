from pydsm.dataset import Dataset
from pydsm import root_sac
import glob
import os

if __name__ == '__main__':
    path = '/work/anselme/TEST_PYDSM/test/201111221848A'
    sac_files = glob.glob(os.path.join(path, '*[RZT]'))

    dataset = Dataset.dataset_from_sac(
        sac_files, verbose=2, headonly=False)
