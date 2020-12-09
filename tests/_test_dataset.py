from dsmpy.dataset import Dataset
from dsmpy import rootdsm_psv
from dsmpy import rootdsm_psv, rootdsm_sh

if __name__ == '__main__':
    parameter_files = [
        rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf']
    dataset = Dataset.dataset_from_files(parameter_files)
    counts, displacements = dataset.get_chunks_station(2)
    counts_eq, displacements_eq = dataset.get_chunks_eq(2)
    
    parameter_files_sh = [rootdsm_sh + 'AK135_SH.inf']
    dataset_sh = Dataset.dataset_from_files(parameter_files_sh, mode=1)
