from dsmpy.seismicmodel import SeismicModel, ModelParameters
from dsmpy.modelparameters import ModelParameters, ParameterType
from dsmpy.dsm import compute_models_parallel
from dsmpy.dataset import Dataset
from dsmpy import root_sac
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys

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

def test_gradient_models():
    model_params = ModelParameters(
        types=[ParameterType.VSH],
        radii=[3480., 3680.],
        mesh_type='boxcar')
    model = SeismicModel.prem().boxcar_mesh(model_params)
    grad_models, dxs = model.gradient_models()

    sac_files = glob.glob(os.path.join(root_sac, '*[RZT]'))
    dataset = Dataset.dataset_from_sac(
        sac_files, headonly=True)

    outputs = compute_models_parallel(
        dataset, grad_models, tlen=1638.4,
        nspc=64, sampling_hz=20, mode=0)
    n_params = len(grad_models) // 2
    n_evs = len(outputs[0])
    waveform_grads = []
    for i in range(n_params):
        waveform_grad_i = []
        for iev in range(len(outputs[0])):
            outputs[i][iev].to_time_domain()
            outputs[i + n_params][iev].to_time_domain()
            waveform_grad_i_iev = (
                outputs[i][iev].us - outputs[i + n_params][iev].us
            ) / (dxs[i] - dxs[i + n_params])
            waveform_grad_i.append(waveform_grad_i_iev)
            outputs[i][iev].free()
            outputs[i + n_params][iev].free()
        waveform_grads.append(waveform_grad_i)
    _, itypes, igrds = model_params.get_free_all_indices()
    types = [model_params.get_types()[i] for i in itypes]
    radii = [model_params.get_grd_params()[i]
             for i in igrds]

    waveform_grads_0000 = np.array([
        5.963111948670274387e-13,7.318878939007591398e-13,
        8.403747805004743159e-13,9.488769106740308336e-13,
        1.057394285868994765e-12,1.165926907905163457e-12,
        1.274474778312821044e-12,1.383037898580892628e-12,
        1.491616270405098464e-12,1.600209895439799779e-12])
    waveform_grads_0010 = np.array([
        -8.391724815037804319e-11, -8.414000083778928765e-11,
        -8.425434896184876587e-11, -8.458559889224371273e-11,
        -8.469999356969719319e-11, -8.492286222369441452e-11,
        -8.492882990515422466e-11, -8.504327156679737037e-11,
        -8.515772889275977600e-11, -8.527220188552298338e-11])
    waveform_grads_0020 = np.array([
        1.229962312554455650e-08, 1.227582010935078955e-08,
        1.224989919567571463e-08, 1.222194127069746893e-08,
        1.219186456734904455e-08, 1.215974998311543457e-08,
        1.212562420531465016e-08, 1.208943257365593747e-08,
        1.205122889372674001e-08, 1.201112123070333273e-08])
    assert np.allclose(waveform_grads_0000, waveform_grads[0][0][0, 0, :10])
    assert np.allclose(waveform_grads_0010, waveform_grads[0][0][1, 0, :10])
    assert np.allclose(waveform_grads_0020, waveform_grads[0][0][2, 0, :10])

if __name__ == '__main__':
    # test_boxcar_mesh()
    test_gradient_models()