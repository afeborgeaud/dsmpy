from obspy import read_events
from pydsm import root_resources
from pydsm.dsm import Event
from pydsm.spc.spctime import SourceTimeFunction
import numpy as np

def convert_catalog():
    cat = read_events(root_resources + 'gcmt.ndk')
    #mts = np.zeros((cat.count(), 6), dtype=np.float64)
    events = np.empty(cat.count(), dtype=np.object)
    for i, event in enumerate(cat):
        tensor = event.preferred_focal_mechanism().moment_tensor.tensor
        stf_obspy = (event.preferred_focal_mechanism().
            moment_tensor.source_time_function)
        source_time_function = SourceTimeFunction(
            stf_obspy.type, 0.5*stf_obspy.duration)
        mt = _mt_from_tensor(tensor)
        lon = event.origins[1].longitude
        lat = event.origins[1].latitude
        depth = event.origins[1].depth / 1000.
        event_id = [e.text for e in event.event_descriptions
                           if e.type == 'earthquake name'][0][1:]
        event = Event(event_id, lat, lon, depth, mt, source_time_function)
        events[i] = event
    np.save(root_resources + 'gcmt', events)

def _mt_from_tensor(tensor):
    mt = np.zeros((3, 3), dtype=np.float64)
    mt[0, 0] = tensor.m_rr
    mt[0, 1] = tensor.m_rt
    mt[0, 2] = tensor.m_rp
    mt[1, 1] = tensor.m_tt
    mt[1, 2] = tensor.m_tp
    mt[2, 2] = tensor.m_pp
    # unit conversion. DSM in units of 10**25 [dyne cm]
    mt *= 1e-18
    return mt
    # return np.array(
    #     [tensor.m_rr, tensor.m_rt, tensor.m_rp,
    #     tensor.m_tt, tensor.m_tp, tensor.m_pp],
    #     dtype=np.float64)

def read_catalog():
    try:
        cat = np.load(root_resources + 'gcmt.npy', allow_pickle=True)
    except:
        convert_catalog()
        cat = np.load(root_resources + 'gcmt.npy', allow_pickle=True)
    return cat

if __name__ == '__main__':
    convert_catalog()