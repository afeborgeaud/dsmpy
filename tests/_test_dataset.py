from pydsm.dataset import Dataset
from pydsm import rootdsm_psv

if __name__ == '__main__':
    parameter_files = [
        rootdsm_psv + 'test2.inf',
        rootdsm_psv + 'test3.inf'
    ]
    dataset = Dataset.dataset_from_files(parameter_files)
    chunks = dataset.get_chunks(4)
    print(chunks)