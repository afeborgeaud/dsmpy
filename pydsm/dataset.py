from pydsm.dsm import PyDSMInput
import numpy as np

class Dataset:
    """Represent a dataset of events and stations.
    """
    def __init__(self, pydsm_inputs):
        self.input_master = pydsm_inputs[0]
        self.lats = np.concatenate([input.lat[:input.nr]
                                    for input in pydsm_inputs])
        self.lons = np.concatenate([input.lon[:input.nr]
                                    for input in pydsm_inputs])
        self.phis = np.concatenate([input.phi[:input.nr]
                                    for input in pydsm_inputs])
        self.thetas = np.concatenate([input.theta[:input.nr]
                                    for input in pydsm_inputs])
        self.eqlats = np.array([input.eqlat for input in pydsm_inputs])
        self.eqlons = np.array([input.eqlon for input in pydsm_inputs])
        self.r0s = np.array([input.r0 for input in pydsm_inputs])
        self.mts = np.concatenate([input.mt for input in pydsm_inputs])
        self.nrs = np.array([input.nr for input in pydsm_inputs])
        self.nr = len(self.lats)

        self.stations = np.concatenate([input.stations
                                        for input in pydsm_inputs])
        self.events = np.array([input.event
                                for input in pydsm_inputs])
        self.source_time_functions = np.array([input.source_time_function
                                              for input in pydsm_inputs])

    @classmethod
    def dataset_from_files(cls, parameter_files):
        pydsm_inputs = [PyDSMInput.input_from_file(file)
                            for file in parameter_files]
        return Dataset(pydsm_inputs)
    

    def get_chunks_station(self, n_cores):
        chunk_size = self.nr // n_cores
        dividers = self.nrs / chunk_size
        dividers = np.round(dividers).astype(np.int64)
        if (dividers == 0).sum() > 0:
            raise RuntimeError('n_cores must be >= number of eqs')
        counts = self.nrs / dividers
        counts_ = []
        for i in range(len(dividers)):
            counts_.append(Dataset._split(self.nrs[i], counts[i]))
        counts = np.concatenate(counts_)
        displacements = counts.cumsum() - counts[0]
        return counts, displacements

    def get_chunks_eq(self, n_cores):
        chunk_size = self.nr // n_cores
        dividers = self.nrs / chunk_size
        dividers = np.round(dividers).astype(np.int64)
        counts = np.ones(dividers.sum())
        displacements_ = []
        for i, divider in enumerate(dividers):
            disp = np.empty(divider)
            disp.fill(i)
            displacements_.append(disp)
        displacements = np.concatenate(displacements_)
        return counts, displacements

    def get_chunks_mt(self, n_cores):
        counts, displacements = self.get_chunks_eq(n_cores)
        return 9*counts, displacements
    
    @staticmethod
    def _split(size, chunk_size):
        n = size // chunk_size
        n = int(n)
        splits = np.empty(n, dtype=np.int64)
        splits.fill(int(chunk_size))
        splits[-1] = size - (n-1) * int(chunk_size)
        return splits


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
