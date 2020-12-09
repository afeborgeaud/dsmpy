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
        self.discon_arr = None
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
            self._n_grd_params = len(radii) * 2 # TODO check consequences
            self._nodes = np.array(radii)
        self.it = 0
        self.set_constraints()
                
    def get_nodes(self):
        return self._nodes

    def get_n_params(self):
        return self._n_grd_params * len(self._types)

    def set_constraints(
            self, mask_dict=None, equal_dict=None,
            discon_arr=None):
        if mask_dict is not None:
            self.mask_dict = mask_dict
        else:
            self.mask_dict = dict()
            for param_type in self._types:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params//2, dtype='bool')
        if equal_dict is not None:
            self.equal_dict = equal_dict
        else:
            self.equal_dict = dict()
            for param_type in self._types:
                self.equal_dict[param_type] = np.arange(
                    self._n_grd_params//2, dtype='int')
        if discon_arr is not None:
            self.discon_arr = discon_arr
        else:
            self.discon_arr = np.ones(
                    self._n_grd_params//2, dtype='bool')

        for param_type in self._types:
            if param_type not in self.equal_dict:
                self.equal_dict[param_type] = np.arange(
                    self._n_grd_params, dtype='int')
            if param_type not in self.mask_dict:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params, dtype='bool')

    def next(self):
        found = False
        istop = 0
        seen_igrd = set()
        while not found or istop < 1000:
            self.it = self.it % self.get_n_params()
            itype = int(self.it // self._n_grd_params)
            igrd = (self.it % self._n_grd_params)
            igrd_2 = igrd // 2
            p_type = self._types[itype]
            if (self.mask_dict[p_type][igrd_2]
                and self.equal_dict[p_type][igrd_2] == igrd_2):
                found = True
            else:
                self.it += 1
            if (igrd_2 not in seen_igrd
                and
                (not self.discon_arr[igrd_2]
                or p_type == ParameterType.RADIUS)):
                self.it += 1
                seen_igrd.add(igrd_2)
                igrd += 1
            istop += 1
        return self.it, itype, igrd

    def get_free_indices(self):
        indices_lst = []
        it_0 = self.it
        self.it = 0
        i = 0
        stop = False
        while not stop:
            i, _, _ = self.next()
            if i in indices_lst:
                stop = True
            else:
                indices_lst.append(i)
            self.it += 1
            i += 1
        self.it = it_0
        return np.array(indices_lst)

    def get_values_matrix(
            self, values_dict):
        values_mat = np.zeros((self._n_grd_params, 9), np.float64)
        for key, values in values_dict.items():
            values_mat[:, key] = values
            
        # TODO check that it's done elsewhere, or modify code
        # TODO check that it doesn't change behavior
        if self.mask_dict is not None:
            for key, mask in self.mask_dict.items():
                mask_expand = np.array(
                    [mask[i//2] for i in range(self._n_grd_params)])
                values_mat[~mask_expand, key] = 0.
        if self.equal_dict is not None:
            for key, indexes in self.equal_dict.items():
                for i, j in enumerate(indexes):
                    if i < j:
                        values_mat[2*i+1, key] = values_mat[2*j+1, key]
        if self.discon_arr is not None:
            for i in np.where(self.discon_arr == False)[0]:
                indices_struct = ParameterType.structure_types()
                values_mat[2*i, indices_struct] = values_mat[
                    2*i+1, indices_struct]
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