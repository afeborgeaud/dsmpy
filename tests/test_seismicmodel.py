from dsmpy.seismicmodel import SeismicModel, ModelParameters
from dsmpy.modelparameters import ModelParameters, ParameterType
import numpy as np
import matplotlib.pyplot as plt

def test_boxcar_mesh():
    model_params = ModelParameters(
        types=[ParameterType.VSH],
        radii=[3480., 3680.],
        mesh_type='boxcar')
    model = SeismicModel.prem().boxcar_mesh(model_params)

    values_dict = {ParameterType.VSH: [0.1]}
    values_1 = model_params.get_values_matrix(values_dict)
    updated_model_1 = model.multiply(values_1)

    values_2 = np.ones(
        (model_params.get_n_grd_params(), 9)) * 0.1
    updated_model_2 = model.multiply(values_2)

    assert np.isclose(
        updated_model_1.get_value_at(3480., ParameterType.VSH),
        7.343071460723877)
    assert np.isclose(
        updated_model_2.get_value_at(3480., ParameterType.VSH),
        7.343071460723877)

if __name__ == '__main__':
    test_boxcar_mesh()
