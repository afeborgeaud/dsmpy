import sys
from pydsm._tish import _tish
from pydsm._tipsv import _pinput, _tipsv
from pydsm.spc import spctime
import numpy as np


class PyDSMOutput:
    def __init__(self, spcs, dsm_input):
        self.spcs = spcs
        self.input = dsm_input
        self.stations = dsm_input.stations
        self.event = dsm_input.event
        self.components = ('Z', 'R', 'T')
        self.dt = 1 / self.input.sampling_hz

    def to_time_domain(self):
        spct = spctime.SpcTime(self.input.tlen, self.input.nspc,
                               self.input.sampling_hz, self.input.omegai,
                               self.input.source_time_function)
        us = spct.spctime(self.spcs)
        self.us = us
        self.ts = np.linspace(0, self.input.tlen,
                              spct.npts, endpoint=False)

    def get_nr(self):
        return self.input.nr


class DSMInput:
    def __init__(
            self, re, ratc, ratl, tlen, nspc, omegai, imin, imax, nzone,
            vrmin, vrmax, rho, vpv, vph, vsv, vsh, eta, qmu, qkappa,
            r0, eqlat, eqlon, mt, nr, theta, phi, lat, lon, output):
        (self.re, self.ratc,
        self.ratl, self.omegai) = re, ratc, ratl, omegai

        (self.tlen, self.nspc,
        self.imin, self.imax) = tlen, nspc, imin, imax

        (self.nzone, self.vrmin, self.vrmax,
        self.rho, self.vpv, self.vph, self.vsv, self.vsh, self.eta,
        self.qmu, self.qkappa) = (
            nzone, vrmin, vrmax, rho, vpv,
            vph, vsv, vsh, eta, qmu, qkappa)

        self.r0, self.eqlat, self.eqlon, self.mt = (
            r0, eqlat, eqlon, mt)

        (self.nr, self.theta, self.phi, self.lat,
        self.lon, self.output) = (nr, theta, phi,
                                 lat, lon, output)

    @classmethod
    def input_from_file(self, parameter_file):
        inputs = _pinput(parameter_file)

        (re, ratc, ratl,
        tlen, nspc, omegai,
        imin, imax, nzone) = inputs[:9]
        (vrmin, vrmax, rho,
        vpv, vph, vsv, vsh,
        eta, qmu, qkappa) = inputs[9:19]
        r0, eqlat, eqlon, mt = inputs[19:23]
        nr = inputs[23]
        (theta, phi, lat,
        lon, output) = inputs[24:]

        return DSMInput(
            re, ratc, ratl, tlen, nspc,
            omegai, imin, imax, nzone, vrmin, vrmax,
            rho, vpv, vph, vsv, vsh, eta, qmu, qkappa, r0, eqlat, eqlon,
            mt, nr, theta, phi, lat, lon, output)

    @classmethod
    def input_from_arrays(self, event, stations,
                          seismicmodel, tlen, nspc):
        pass  # TODO

    def get_inputs_for_tish(self):
        # TODO modify fortran? Else, have to take care of case
        # number of core layers != 2
        nzone = self.nzone - 2
        vrmin = np.pad(self.vrmin[2:], (0,2), constant_values=0)
        vrmax = np.pad(self.vrmax[2:], (0,2), constant_values=0)
        qmu = np.pad(self.qmu[2:], (0,2), constant_values=0)
        npad = ((0, 0), (0, 2))
        rho = np.pad(self.rho[:,2:], npad, constant_values=0)
        vsv = np.pad(self.vsv[:,2:], npad, constant_values=0)
        vsh = np.pad(self.vsh[:,2:], npad, constant_values=0)

        inputs = (self.re, self.ratc, self.ratl, self.tlen,
                  self.nspc, self.omegai, self.imin, self.imax,
                  nzone, vrmin, vrmax, rho,
                  vsv, vsh, qmu, self.r0,
                  self.eqlat, self.eqlon, self.mt, self.nr,
                  self.theta, self.phi, self.lat, self.lon, self.output)
        return inputs

    def get_inputs_for_tipsv(self):
        inputs = (
            self.re, self.ratc, self.ratl, self.tlen,
            self.nspc, self.omegai, self.imin, self.imax,
            self.nzone, self.vrmin, self.vrmax, self.rho,
            self.vpv, self.vph, self.vsv, self.vsh, self.eta,
            self.qmu, self.qkappa, self.r0, self.eqlat, self.eqlon,
            self.mt, self.nr, self.theta, self.phi, self.lat, self.lon,
            self.output)
        return inputs


class PyDSMInput(DSMInput):
    def __init__(
            self, dsm_input, sampling_hz=None,
            source_time_function=None):
        super().__init__(
            *dsm_input.get_inputs_for_tipsv())
        self.source_time_function = source_time_function
        self.sampling_hz = self.find_optimal_sampling_hz(sampling_hz)
        self.stations = self._parse_stations()
        self.event = self._parse_event()

    @classmethod
    def input_from_file(self, parameter_file,
                        sampling_hz=None, source_time_function=None):
        dsm_input = super().input_from_file(parameter_file)
        pydsm_input = PyDSMInput(dsm_input, sampling_hz,
                                 source_time_function)
        return pydsm_input

    @classmethod
    def input_from_arrays(
            self, event, stations,
            seismicmodel, tlen, nspc, source_time_function,
            sampling_hz):
        dsm_input = super().input_from_arrays(event, stations,
                                              seismicmodel, tlen, nspc)
        pydsm_input = PyDSMInput(dsm_input, sampling_hz,
                                 source_time_function)
        return pydsm_input

    def find_optimal_sampling_hz(self, sampling_hz):
        if sampling_hz is not None:
            return sampling_hz
        else:
            Tmin = (self.tlen
                    / self.imax)
            optimal_sampling_hz = 40 / Tmin  # TODO: check this scaling
            return optimal_sampling_hz

    def set_sourcetimefunction(self, sourcetimefunction):
        self.sourcetimefunction = sourcetimefunction

    def _parse_stations(self):
        stations = []
        for i in range(self.nr):
            name, net = self.output[i].tostring(). \
                decode('utf-8').split('/')[-1].split('.')[0].split('_')
            station = Station(name, net, self.lat[i], self.lon[i])
            stations.append(station)
        return tuple(stations)

    def _parse_event(self):
        eventID = self.output[0].tostring(). \
            decode('utf-8').split('/')[-1].split('.')[1]
        if eventID[-2:] == 'SH':
            eventID = eventID[:-2]
        elif eventID[-3:] == 'PSV':
            eventID = eventID[:-3]
        else:
            raise RuntimeError('{}'.format(eventID))
        event = Event(eventID, self.eqlat, self.eqlon,
                      6371. - self.r0, self.mt)
        return event


class Station:
    def __init__(self, name: str, network: str,
            latitude: float, longitude: float):
        self.name = name
        self.network = network
        self.latitude = latitude
        self.longitude = longitude

    def __repr__(self):
        return self.name + '_' + self.network


class Event:
    """Represent an earthquake point-source.

    Attributes
    ----------
    eventID: str
        GCMT event name
    latitude: float
        centroid geographic latitude [-90, 90] in degree
    longitude: float
        centroid longitude [-180, 180] in degree
    depth: float
        centroid depth in km
    mt: ndarray(3, 3)
        moment tensor
    """

    def __init__(self, eventID, latitude, longitude, depth, mt):
        self.eventID = eventID
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.mt = mt


def compute(dsm_input, mode=0, write_to_file=False):
    """Compute spectra using DSM.

    Args:
        dsm_input (DSMInput): inputs for DSM
        mode (int): computation mode. 0: both, 1: P-SV, 2: SH
        write_to_file (bool): write spetrum files to disk
            as specified in dsm_input.output

    Returns:
        dsm_output (DSMOutput): object containing spectra and
            statations/source information
    
    Note:
        SH and P-SV spectra are summed by default. Using only P-SV
        or SH results in non-physical waves and should be avoided.
        See Kawai et al. (2006) for details.
    """
    print('compute SH')
    sh_spcs = _tish(*dsm_input.get_inputs_for_tish(),
        write_to_file)
    # FIXME memory error in tipsv.tipsv
    #print('compute PSV')
    #psv_spcs = _tipsv(*dsm_input.get_inputs_for_tipsv(),
    #    write_to_file)
    #spcs = sh_spcs + psv_spcs
    spcs = sh_spcs
    dsm_output = PyDSMOutput(spcs, dsm_input)
    return dsm_output
