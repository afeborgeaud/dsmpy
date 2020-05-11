'''Interface for Fortran tish
'''

from pydsm.lib import tish


def _pinput(parameter_file):
    inputs = tish.pinput_fromfile(parameter_file)
    return inputs


def _tish(
        re, ratc, ratl, tlen, nspc, omegai,
        imin, imax, nzone, vrmin, vrmax, rho,
        vsv, vsh, qmu, r0, eqlat, eqlon, mt, nr,
        theta, phi, lat, lon, output, write_to_file):
    spcs = tish.tish(
        re, ratc, ratl, tlen, nspc, omegai,
        imin, imax, nzone, vrmin, vrmax, rho,
        vsv, vsh, qmu, r0, eqlat, eqlon, mt, nr,
        theta, phi, lat, lon, output, write_to_file)
    return spcs
