from pydsm.dsm import PyDSMInput
import numpy as np

class Dataset:
    """Represent a dataset of events and stations.
    """

    def __init__(self, pydsm_inputs):
        self.pydsm_inputs = pydsm_inputs
        self.nr = np.sum([len(input.stations) 
                          for input in pydsm_inputs])
        self.ne = len(pydsm_inputs)

    @classmethod
    def dataset_from_files(cls, parameter_files):
        pydsm_inputs = [PyDSMInput.input_from_file(file)
                            for file in parameter_files]
        return Dataset(pydsm_inputs)
    
    # def split(self, chunk_size=100):
    #     inputs = []
    #     for input in self.pydsm_inputs:
    #         nr_ = input.nr // chunk_size
    #         if n_new % chunk_size > 0:
    #             nr_ += 1
    #         for i in range(nr_ - 1):
                
    #             input_ = PyDSMInput(input.dsm_input, )