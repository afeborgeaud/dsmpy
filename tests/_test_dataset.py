from pydsm.dataset import Dataset
from pydsm import rootdsm_psv

if __name__ == '__main__':
    from pydsm import rootdsm_psv
    parameter_files = [
        rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf']
    dataset = Dataset.dataset_from_files(parameter_files)
    counts, displacements = dataset.get_chunks_station(2)
    counts_eq, displacements_eq = dataset.get_chunks_eq(2)
    print(counts)
    print(displacements.dtype)