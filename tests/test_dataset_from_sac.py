from pydsm.dataset import Dataset
from pydsm import root_sac
import glob
import os

if __name__ == '__main__':
    sac_files = glob.glob(os.path.join(root_sac, '*T'))
    
    dataset_info = Dataset.dataset_from_sac(sac_files)
