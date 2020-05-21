from pydsm.dataset import Dataset
from pydsm import root_sac
import glob
import os

if __name__ == '__main__':
    path = '/work/anselme/TEST_PYDSM/test/200609220232A'
    sac_files = glob.glob(os.path.join(root_sac, '*T'))
    
    dataset = Dataset.dataset_from_sac(sac_files)
