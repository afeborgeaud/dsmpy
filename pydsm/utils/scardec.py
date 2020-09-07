from pydsm import root_resources
import numpy as np
import glob
from pydsm.utils.cmtcatalog import read_catalog
from pydsm.event import Event
from pydsm.spc.stf import SourceTimeFunction
from pydsm.spc.stfcatalog import STFCatalog
import matplotlib.pyplot as plt
import os

def get_stf(event):
    '''Return source time function in time domain.
    Args:
        event (pydsm.Event): event
    Returns:
        stf (ndarray((2,npts))): source time function, normalized so
            that its integral is 1
    '''
    dir_stf = _parse_dir_name(event)
    if dir_stf is None:
        return None

    file_name = glob.glob(dir_stf + '/fctoptsource*')[0]
    stf = np.loadtxt(file_name, skiprows=2)
    start, end = _get_start_end(stf)
    stf = stf[start:end]
    stf_integral = _compute_integral(stf)
    stf[:, 1] /= stf_integral
    stf[:, 0] -= stf[0, 0]
    return stf

def get_duration(event):
    stf = get_stf(event)
    if stf is None:
        return None
    return stf[-1, 0]

def create_catalog():
    cmt_catalog = read_catalog()
    stf_catalog = dict()
    for event in cmt_catalog:
        print(event)
        dir_name = _parse_dir_name(event)
        if dir_name:
            duration = get_duration(event)
            stf = SourceTimeFunction(
                'triangle', duration/2.)
            stf_catalog[event.event_id] = stf

    path = os.path.joint(root_resources, 'scardec.pkl')
    STFCatalog.save(path, scardec_catalog)

def _parse_dir_name(event):
    dir_scardec = os.path.join(root_resources, 'scardec')
    # event_id_post2005 = _convert_name_to_post2005(event)
    partial_scardec_dir_name = _convert_name_to_partial_scardec(event)
    dirs = glob.glob(dir_scardec + '/*' + partial_scardec_dir_name + '*')
    print(dirs)
    parsed_dir = None
    for dir_ in dirs:
        ss = int(dir_.split(partial_scardec_dir_name)[1][:2])
        if np.abs(ss - event.centroid_time.second) <= 15:
            parsed_dir = dir_
    return parsed_dir

def _convert_name_to_partial_scardec(event):
    if event.centroid_time is not None:
        return event.centroid_time.strftime('%Y%m%d_%H%M')

def _convert_name_to_post2005(event):
    if event.centroid_time is not None:
        return event.centroid_time.strftime('%Y%m%d%H%M')
    else:
        return None

def _get_start_end(stf):
    i_peak = np.argmax(stf[:, 1])

    i = i_peak
    while (stf[i, 1] > 0) and (i < len(stf)):
        i += 1
    i_end = i + 1

    i = i_peak
    while (stf[i, 1] > 0) and (i >= 0):
        i -= 1
    i_start = i

    return i_start, i_end
    

def _compute_integral(stf):
    integral = np.trapz(stf[:,1], stf[:,0])
    return integral

if __name__ == '__main__':
    # create_catalog()
    