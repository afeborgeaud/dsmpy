import sys
from pyDSM.lib import tish
from pyDSM.spc import spctime
import numpy as np

class DSMoutput:
    def __init__(self, spcs, dsm_input):
        self.spcs = spcs
        self.input = dsm_input
        self.stations = dsm_input.input_parameters.get_station_list()
        self.event = dsm_input.input_parameters.get_event()
    
    def to_time_domain(self):
        spct = spctime.SpcTime(self.input.get_tlen(), 
            self.input.get_nspc(), self.input.samplingHz, 
            self.input.get_omegai(), sourcetimefunction=self.input.sourcetimefunction)
        u = spct.spctime(self.spcs)
        self.u = u
        self.ts = np.linspace(0, self.input.get_tlen(),
            spct.npts, endpoint=False)

    def get_nr(self):
        return self.input.input_parameters.nr

class InputParameters:
    def __init__(self, inputs):
        self.re, self.ratc, self.ratl, \
        self.tlen, self.np, self.omegai, \
        self.imin, self.imax, self.nzone = inputs[:9]

        self.vrmin, self.vrmax, self.rho, \
        self.vsv, self.vsh, self.qmu = inputs[9:15]
        #[x[:self.nzone] for x in inputs[9:15]]

        self.r0, self.eqlat, self.eqlon, self.mt = inputs[15:19]

        self.nr = inputs[19]
        self.theta, self.phi, self.lat, \
        self.lon, self.output = inputs[20:]
        #[x[:self.nr] for x in inputs[20:]]
    
    def get_inputs(self):
        inputs = (self.re, self.ratc, self.ratl, self.tlen,
            self.np, self.omegai, self.imin, self.imax,
            self.nzone, self.vrmin, self.vrmax, self.rho,
            self.vsv, self.vsh, self.qmu, self.r0,
            self.eqlat, self.eqlon, self.mt, self.nr,
            self.theta, self.phi, self.lat, self.lon, self.output)
        return inputs

    def get_station_list(self):
        stations = []
        for i in range(self.nr):
            name, net = self.output[i].tostring().\
                decode('utf-8').split('/')[-1].split('.')[0].split('_')
            station = Station(name, net, self.lat[i], self.lon[i])
            stations.append(station)
        return tuple(stations)

    def get_event(self):
        eventID = self.output[0].tostring().\
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

class DSMinput:
    def __init__(self, parameterfile, samplingHz=None):
        self.parameterfile = parameterfile
        self.sourcetimefunction = None
        self.input_parameters = self.pinput()
        self.set_samplingHz(samplingHz)

    def set_samplingHz(self, samplingHz):
        if samplingHz is not None:
            self.samplingHz = samplingHz
        else:
            Tmin = self.input_parameters.tlen / self.input_parameters.imax
            self.samplingHz = 40 / Tmin # TODO: check this scaling

    def set_sourcetimefunction(self, sourcetimefunction):
        self.sourcetimefunction = sourcetimefunction
    
    def pinput(self):
        inputs = tish.pinput_fromfile(self.parameterfile)
        return InputParameters(inputs)

    def get_tlen(self):
        return self.input_parameters.tlen

    def get_nspc(self):
        return self.input_parameters.np
    
    def get_omegai(self):
        return self.input_parameters.omegai

class Station:
    def __init__(self, name, network, latitude, longitude):
        self.name = name
        self.network = network
        self.latitude = latitude
        self.longitude = longitude
    def __repr__(self):
        return self.name + '_' + self.network
    
class Event:
    def __init__(self, eventID, latitude, longitude, depth, mt):
        self.eventID = eventID
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.mt = mt

def compute(dsm_input, write_to_file=False):
    spcs = tish.tish(*dsm_input.input_parameters.get_inputs(),
        write_to_file)
    dsm_output = DSMoutput(spcs, dsm_input)
    return dsm_output