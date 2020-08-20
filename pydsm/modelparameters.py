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
        if mesh_type == 'boxcar':
            self._n_nodes = len(radii) - 1
            self._nodes = np.array(radii)
        elif mesh_type == 'triangle':
            # self._n_nodes = int(
            #     len(radii) - 1 + np.ceil((len(radii)-1) / 2))
            self._n_nodes = 2 * len(radii) - 1
            n = 2 * len(radii) - 1
            self._nodes = np.zeros(n, dtype=np.float64)
            for i in range(n-1):
                self._nodes[i] = (radii[i//2] + (i % 2)
                    * (radii[i//2+1] - radii[i//2]) / 2.)
            self._nodes[-1] = radii[-1]
                
    def get_nodes(self):
        return self._nodes

    def get_values_matrix(self, values_dict):
        values_mat = np.zeros((self._n_nodes, 8), np.float64)
        for key, values in values_dict.items():
            values_mat[:, key] = values
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

if __name__ == '__main__':
    types = [ParameterType.VSV, ParameterType.VSH]
    radii = np.array([3480., 3700.], dtype=np.float64)
    model_params = ModelParameters(types, radii)
    values_dict = {
        ParameterType.VSV: -1.,
        ParameterType.VSH: 1.}
    values_mat = model_params.get_values_matrix(values_dict)
    print(values_mat)