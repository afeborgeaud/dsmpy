import numpy as np
from enum import IntEnum

class ModelParameters:
    """Represent model parameters.

    Args:
        types (list(modelparameters.ParameterType)): e.g., RHO, VSH
        radii (ndarray(float)): radii of boundaries of perturbed layers
        mesh_type (str): 'boxcar' or 'triangle'
    """
    def __init__(self, types, radii, mesh_type='boxcar'):
        self._types = types
        self._radii = radii
        self._mesh_type = mesh_type
        self.mask_dict = None
        self.equal_dict = None
        if mesh_type == 'boxcar':
            self._n_nodes = len(radii)
            self._n_grd_params = len(radii) - 1
            self._nodes = np.array(radii)
        elif mesh_type == 'triangle':
            # self._n_nodes = int(
            #     len(radii) - 1 + np.ceil((len(radii)-1) / 2))
            self._n_nodes = len(radii)
            self._n_grd_params = 2 * len(radii) - 1
            n = 2 * len(radii) - 1
            self._nodes = np.zeros(n, dtype=np.float64)
            for i in range(n-1):
                self._nodes[i] = (radii[i//2] + (i % 2)
                    * (radii[i//2+1] - radii[i//2]) / 2.)
            self._nodes[-1] = radii[-1]
        elif mesh_type == 'lininterp':
            self._n_nodes = len(radii)
            self._n_grd_params = len(radii)
            self._nodes = np.array(radii)
        self.it = 0
                
    def get_nodes(self):
        return self._nodes

    def get_n_params(self):
        return self._n_grd_params * len(self._types)

    def set_constraints(self, mask_dict=None, equal_dict=None):
        if mask_dict is not None:
            self.mask_dict = mask_dict
        else:
            self.mask_dict = dict()
            for param_type in self._types:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params, dtype='bool')
        if equal_dict is not None:
            self.equal_dict = equal_dict
        else:
            self.equal_dict = dict()
            for param_type in types:
                self.equal_dict[param_type] = np.arange(
                    model_params._n_grd_params, dtype='int')

        for param_type in self._types:
            if param_type not in self.equal_dict:
                self.equal_dict[param_type] = np.arange(
                    self._n_grd_params, dtype='int')
            if param_type not in self.mask_dict:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params, dtype='bool')

    def get_it_indices(self):
        found = False
        istop = 0
        while not found or istop == 1000:
            self.it = self.it % (self._n_grd_params*len(self._types))
            itype = int(self.it // self._n_grd_params)
            igrd = self.it % self._n_grd_params
            if (self.mask_dict[self._types[itype]][igrd]
                and self.equal_dict[self._types[itype]][igrd] == igrd):
                found = True
            else:
                self.it += 1
            istop += 1
        return self.it, itype, igrd

    def get_free_indices(self):
        mask_arr = np.ones(self.get_n_params(), dtype='bool')
        equal_arr = np.array(mask_arr)
        if self.mask_dict is not None:
            mask_arr = np.hstack(
                [self.mask_dict[p] for p in self._types])
            print(mask_arr)
        if self.equal_dict is not None:
            equal_arr = np.hstack(
                [(self.equal_dict[p]
                ==np.arange(len(self.equal_dict[p]), dtype='int')) 
                for p in self._types])
            print(equal_arr)
        indices = np.where(mask_arr & equal_arr)[0]
        return indices

    def get_values_matrix(
            self, values_dict):
        values_mat = np.zeros((self._n_grd_params, 9), np.float64)
        for key, values in values_dict.items():
            values_mat[:, key] = values
            
        # TODO check that it's done elsewhere, or modify code structure
        # if self.mask_dict is not None:
        #     for key, mask in self.mask_dict.items():
        #         values_mat[~mask, key] = 0.
        # if self.equal_dict is not None:
        #     for key, indexes in self.equal_dict.items():
        #         for i, j in enumerate(indexes):
        #             values_mat[i, key] = values_mat[j, key]
        return values_mat


class ParameterType(IntEnum):
    RHO = 0
    VPV = 1
    VPH = 2
    VSV = 3
    VSH = 4
    ETA = 5
    QMU = 6
    QKAPPA = 7
    RADIUS = 8

    @staticmethod
    def structure_types():
        return [ParameterType(i) for i in range(6)]

if __name__ == '__main__':
    types = [ParameterType.VSV, ParameterType.VSH]
    radii = np.array([3480., 3700.], dtype=np.float64)
    model_params = ModelParameters(types, radii)
    values_dict = {
        ParameterType.VSV: -1.,
        ParameterType.VSH: 1.}
    values_mat = model_params.get_values_matrix(values_dict)
    print(values_mat)