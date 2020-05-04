import sys
from pyDSM.lib import tish
from pyDSM.spc import spctime
import numpy as np

class SourceTimeFunction:
    @staticmethod
    def triangle(nspc):
        stf = np.array((nspc), dtype=np.complex128)
        return stf

class DSMoutput:
    def __init__(self, spcs, dsm_input):
        self.spcs = spcs
        self.input = dsm_input
    
    def to_time_domain(self):
        spct = spctime.SpcTime(self.input.get_tlen(), 
            self.input.get_nspc(), self.input.samplingHz, 
            self.input.get_omegai())
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

class DSMinput:
    def __init__(self, parameterfile, sourcetimefunction=None, samplingHz=20):
        self.parameterfile = parameterfile
        self.sourcetimefunction = sourcetimefunction
        self.samplingHz = samplingHz
        self.input_parameters = self.pinput()
    
    def pinput(self):
        inputs = tish.pinput_fromfile(self.parameterfile)
        return InputParameters(inputs)

    def get_tlen(self):
        return self.input_parameters.tlen

    def get_nspc(self):
        return self.input_parameters.np
    
    def get_omegai(self):
        return self.input_parameters.omegai
    

def compute(dsm_input, write_to_file=False):
    spcs = tish.tish(*dsm_input.input_parameters.get_inputs(),
        write_to_file)
    dsm_output = DSMoutput(spcs, dsm_input)
    return dsm_output