'''Interface for Fortran tish
'''

from pydsm.lib import tish
import numpy as np

parameters = dict(maxnzone=tish.parameters.maxnzone.item(),
                  maxnr=tish.parameters.maxnr.item(),
                  maxlmax=tish.parameters.maxlmax.item(),
                  maxnlay=tish.parameters.maxnlay.item(),
                  flattening=tish.parameters.flattening.item())

def _pinput(parameter_file):
    """DSM input from DSM input parameter file.

    Args:
        parameter_file (str): path to DSM parameter file
            (e.g. 'PREM_SH.inf')
    Returns:
        inputs (tuple): tupe of input parameters 
            (in correct order for _tish)
    """
    inputs = tish.pinput_tish(parameter_file)
    return inputs

def _tish(
        re, ratc, ratl, tlen, nspc, omegai,
        imin, imax, nzone, vrmin, vrmax, rho,
        vsv, vsh, qmu, r0, eqlat, eqlon, mt, nr,
        theta, phi, lat, lon, output, write_to_file):
    """Compute spectra for SH wavefield using DSM."""
    spcs = tish.tish(
        re, ratc, ratl, tlen, nspc, omegai,
        imin, imax, nzone, vrmin, vrmax, rho,
        vsv, vsh, qmu, r0, eqlat, eqlon, mt, nr,
        theta, phi, lat, lon, output, write_to_file)
    return spcs

def _translat(lat_geodetic):
    """Return geocentric latitude from geodetic latitude.
    Notes:
        flattening = 1/298.25
    """
    return tish.translat(lat_geodetic)

def _calthetaphi(stalat, stalon, eqlat, eqlon):
    """Compute theta and phi for spherical coordinate system.

    Args:
        eqlat (float): event geodetic latitude
        eqlon (float): event longitude
        stalat (float): station geodetic latitude
        stalon (float): station longitude
    Returns:
        theta (float): theta spherical coordinate
        phi (float): phi spherical coordinate
    """
    eqlat_geocentric = _translat(eqlat)
    stalat_geocentric = _translat(stalat)
    theta, phi = tish.calthetaphi(eqlat_geocentric, eqlon,
                            stalat_geocentric, stalon)
    # the conversion to radians is now done in pinput
    # so it must be reproduced here
    theta *= np.pi / 180.
    phi *= np.pi / 180.
    return theta, phi

