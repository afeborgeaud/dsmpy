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
    model_ref = ak135_lin_dpp()
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

def ak135_lin_dpp():
    model = SeismicModel.ak135()

    r0 = 4000.
    model = model._add_boundary(r0)
    izone = model.get_zone(r0 - 1e-5)
    vp_cmb = model.get_value_at(3479.5, ParameterType.VPH)
    vp_r0 = model.get_value_at(r0, ParameterType.VPH)
    vp_cmb_p = vp_cmb + 0.09
    vs_cmb = model.get_value_at(3479.5, ParameterType.VSH)
    vs_r0 = model.get_value_at(r0, ParameterType.VSH)

    x_cmb = 3479.5 / 6371
    x0 = r0 / 6371
    poly_vp = model._lin_element(x_cmb, x0,
                                     vp_cmb_p, vp_r0)
    poly_vs = model._lin_element(x_cmb, x0,
                                     vs_cmb, vs_r0)

    for p_type in [ParameterType.VPH, ParameterType.VPV]:
        model.set_value(izone, p_type, poly_vp)
        model.set_value(izone - 1, p_type, poly_vp)
    for p_type in [ParameterType.VSH, ParameterType.VSV]:
        model.set_value(izone, p_type, poly_vs)
        model.set_value(izone - 1, p_type, poly_vs)

    return model


def get_range_dict(model_params, delta_dict):
    range_dict = dict()
    for param_type in model_params._types:
        range_arr = np.empty((model_params._n_grd_params, 2), dtype='float')
        delta = delta_dict[param_type]
        range_arr[:, 0] = -delta
        range_arr[:, 1] = delta
        range_dict[param_type] = range_arr

    return range_dict

if __name__ == '__main__':
    single_layer_dpp()