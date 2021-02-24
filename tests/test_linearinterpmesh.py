"""Test script for building seismic models using
a discontinuous linear spline parametrization.
"""

from dsmpy.modelparameters import ModelParameters, ParameterType
from dsmpy.seismicmodel import SeismicModel
from dsmpy.utils.modelutils import ak135_lin_dpp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    types = [ParameterType.VSH, ParameterType.RADIUS]
    radii = 6371. - np.array([2891.5, 2691.5])
    model_ref = ak135_lin_dpp()

    model_params = ModelParameters(
        types, radii, mesh_type='lininterp')
    model = model_ref.lininterp_mesh(
        model_params, discontinuous=False)

    # constraints to parameters
    mask_dict = dict()
    mask_dict[ParameterType.VSH] = np.ones(
        model_params._n_grd_params, dtype='bool')
    mask_dict[ParameterType.RADIUS] = np.ones(
        model_params._n_grd_params, dtype='bool')
    mask_dict[ParameterType.RADIUS][[0, 1]] = False

    equal_dict = dict()
    equal_dict[ParameterType.VSH] = np.arange(
        model_params._n_grd_params, dtype='int')
    equal_dict[ParameterType.VSH][2] = 0
    equal_dict[ParameterType.VSH][3] = 1

    discon_arr = np.zeros(
        model_params._n_nodes, dtype='bool')
    discon_arr[1] = True

    model_params.set_constraints(
        mask_dict=mask_dict,
        equal_dict=equal_dict,
        discon_arr=discon_arr)

    print('Free indices:', model_params.get_free_indices())

    rng = np.random.default_rng(seed=0)
    fig, axes = plt.subplots(1, 3, figsize=(12, 6))
    for i in range(10):
        values = np.zeros((model_params._n_grd_params, len(types)))
        if i > 0:
            values[:, 0] = rng.uniform(-0.5, 0.5, size=(values.shape[0]))
            values[:, 1] = rng.uniform(-190, 190, size=(values.shape[0]))
        values_dict = {
            param_type: values[:, i] for i, param_type in enumerate(types)
        }
        model_gen = model.build_model(model, model_params, values_dict)

        model_gen.plot(ax=axes[0], types=[ParameterType.VSH])
        model_gen.plot(ax=axes[1], types=[ParameterType.VPH])
        model_gen.plot(ax=axes[2], types=[ParameterType.RHO])

    axes[0].set_xlim([6.5, 8.])
    axes[1].set_xlim([13, 14.])
    axes[2].set_xlim([5.2, 5.8])
    for ax in axes:
        ax.set_ylim([3480., 4200.])

    fig.savefig('./figures/models_lininterp.pdf', bbox_inches='tight')