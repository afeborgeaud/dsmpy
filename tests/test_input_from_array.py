from dsmpy import rootdsm_sh
from dsmpy._tish import _tish, _pinput
from dsmpy.dsm import DSMInput, PyDSMInput
from dsmpy.station import Station
from dsmpy.event import Event, MomentTensor
from dsmpy.seismicmodel import SeismicModel
from dsmpy.dataset import Dataset
import numpy as np
import os
from datetime import datetime
import joblib
import glob


def _get_ref_pydsm_input():
    parameter_file = os.path.join(rootdsm_sh, 'AK135_SH_64.inf')
    inputs = PyDSMInput.input_from_file(parameter_file, file_mode=2)
    return inputs


def test_input_from_arrays():
    seismic_model = SeismicModel.ak135()

    mt = MomentTensor(Mrr=-0.193, Mrt=0.0998, Mrp=-0.148,
                                 Mtt=0.0571, Mtp=-0.141, Mpp=0.136)
    source_time_function = None

    event = Event('', -1.600000023841858, -78.05000305175781,
                      6371-6195.2, mt,
                      centroid_time=datetime(1999, 1, 1),
                      source_time_function=source_time_function)
    stations = [
        Station('109C', 'TA', 32.88890075683594, -117.1051025390625)]
    tlen = 3276.8
    nspc = 64
    sampling_hz = 20

    pydsm_input = PyDSMInput.input_from_arrays(
        event, stations, seismic_model,
        tlen, nspc, sampling_hz)
    
    # reference input to test against
    ref_input = _get_ref_pydsm_input()

    # input from array return symetric mt
    # TODO check why DSM works with pinput returning
    # non-symetric mt
    mt_ = ref_input.mt
    mt_[1,0] = mt_[0,1]
    mt_[2,0] = mt_[0,2]
    mt_[2,1] = mt_[1,2]

    assert pydsm_input.omegai == ref_input.omegai
    assert pydsm_input.re == ref_input.re
    assert pydsm_input.ratc == ref_input.ratc
    assert pydsm_input.ratl == ref_input.ratl
    assert pydsm_input.tlen == ref_input.tlen
    assert pydsm_input.nspc == ref_input.nspc
    assert pydsm_input.imin == ref_input.imin
    assert pydsm_input.imax == ref_input.imax
    assert pydsm_input.lat[0] == ref_input.lat[0]
    assert pydsm_input.lon[0] == ref_input.lon[0]
    assert pydsm_input.theta[0] == ref_input.theta[0]
    assert pydsm_input.phi[0] == ref_input.phi[0]
    assert np.allclose(pydsm_input.mt, mt_, atol=1e-10)
    assert np.allclose(pydsm_input.get_inputs_for_tish()[11], ref_input.rho)
    print('All passed!')


def test_input_for_dataset_parallel():
    sac_files = list(glob.iglob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sac_files_2/*')))[:3]

    # Create the dataset
    dataset = Dataset.dataset_from_sac(sac_files, headonly=False)

    dataset.stations

    # define computation parameters
    tlen = 1638.4
    nspc = 64
    sampling_hz = 20
    mode = 2

    model = SeismicModel.prem()

    input = PyDSMInput.input_from_arrays(
        dataset.events[0], dataset.stations, model,
        tlen, nspc, sampling_hz)

    in_tish = input.get_inputs_for_tish()

    parameter_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'input_files/AK135_SH_64.inf')
    ref_inputs = _pinput(parameter_file)

    assert len(in_tish) == len(ref_inputs)
    for i in range(len(in_tish)):
        if type(in_tish[i]) in [int, float, np.float64]:
            pass
        else:
            assert in_tish[i].shape == ref_inputs[i].shape

    print('Solving with ref_inputs')
    _tish(*ref_inputs, write_to_file=False)
    print('Solving with inputs from arrays')
    _tish(*in_tish, write_to_file=False)


if __name__ == '__main__':
    test_input_from_arrays()
    test_input_for_dataset_parallel()
    