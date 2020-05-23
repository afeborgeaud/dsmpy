import numpy as np
from enum import IntEnum

class ModelParameters:
    """Represent model parameters.
    """
    def __init__(self, types, radii):
        self._types = types
        self._radii = radii
        self._n_nodes = len(radii)
    
    def get_nodes(self):
        return self._radii

    def get_values_matrix(self, values_dict):
        values_mat = np.zeros((self._n_nodes-1, 8))
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