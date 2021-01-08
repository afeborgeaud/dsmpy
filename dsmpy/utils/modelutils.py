"""Utilities to build various model meshes."""

from dsmpy.seismicmodel import SeismicModel
from dsmpy.modelparameters import ModelParameters, ParameterType
import numpy as np
import matplotlib.pyplot as plt

def single_layer_dpp():
    """Create objects for a single-layer D'' model.

    Returns:
        SeismicModel: reference model
        ModelParameters: model parameters
        dict: range dict

    """
    # model parameters
    model_ref = SeismicModel.ak135()
    types = [ParameterType.VSH, ParameterType.RADIUS]
    radii = np.array([3479.5, 3679.5])
    model_params = ModelParameters(
        types, radii, mesh_type='lininterp')

    # Defines constraints to model parameters
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

    # Defines value range for model parameters
    delta_dict = {
        ParameterType.VSH: 0.5,
        ParameterType.RADIUS: 190.}
    range_dict = get_range_dict(model_params, delta_dict)

    return model_ref, model_params, range_dict


def get_range_dict(model_params, delta_dict):
    range_dict = dict()
    for param_type in model_params._types:
        range_arr = np.empty((model_params._n_grd_params, 2), dtype='float')
        delta = delta_dict[param_type]
        range_arr[:, 0] = -delta
        range_arr[:, 1] = delta
        range_dict[param_type] = range_arr

    return range_dict

single_layer_dpp()