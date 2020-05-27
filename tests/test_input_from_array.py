from pydsm import dsm, rootdsm_sh
from pydsm.station import Station
from pydsm.event import Event, MomentTensor
from pydsm.seismicmodel import SeismicModel
import numpy as np
import os

def get_ref_pydsm_input():
    parameter_file = os.path.join(rootdsm_sh, 'AK135_SH_64.inf')
    inputs = dsm.PyDSMInput.input_from_file(parameter_file, mode=2)
    return inputs

if __name__ == '__main__':
    seismic_model = SeismicModel.ak135()

    mt = MomentTensor(Mrr=-0.193, Mrt=0.0998, Mrp=-0.148,
                                 Mtt=0.0571, Mtp=-0.141, Mpp=0.136)
    source_time_function = None

    event = Event('', -1.600000023841858, -78.05000305175781,
                      6371-6195.2, mt.to_array(),
                      source_time_function=source_time_function)
    stations = [
        Station('109C', 'TA', 32.88890075683594, -117.1051025390625)]
    tlen = 3276.8
    nspc = 64
    sampling_hz = 20

    pydsm_input = dsm.PyDSMInput.input_from_arrays(
        event, stations, seismic_model,
        tlen, nspc, sampling_hz)
    
    # reference input to test against
    ref_input = get_ref_pydsm_input()

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
    