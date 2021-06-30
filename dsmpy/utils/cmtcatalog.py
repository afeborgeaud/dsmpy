from dsmpy import root_resources
from dsmpy.event import Event, MomentTensor
from dsmpy.spc.stf import SourceTimeFunction
from obspy import read_events
import numpy as np
import warnings
from datetime import date
import re
import requests

def convert_catalog(cat):
    #cat = read_events(root_resources + 'gcmt.ndk')
    #mts = np.zeros((cat.count(), 6), dtype=np.float32)
    events = np.empty(cat.count(), dtype=np.object)
    for i, event in enumerate(cat):
        tensor = event.preferred_focal_mechanism().moment_tensor.tensor
        mw = [m for m in event.magnitudes 
              if m.magnitude_type=='Mwc'][0].mag
        centroid_time = [o for o in event.origins
                        if o.origin_type == 'centroid'][0].time
        stf_obspy = (event.preferred_focal_mechanism().
            moment_tensor.source_time_function)
        source_time_function = SourceTimeFunction(
            stf_obspy.type, 0.5*stf_obspy.duration)
        mt = _mt_from_tensor(tensor, mw)
        lon = event.origins[1].longitude
        lat = event.origins[1].latitude
        depth = event.origins[1].depth / 1000.
        event_id = [e.text for e in event.event_descriptions
                           if e.type == 'earthquake name'][0][1:]
        event = Event(event_id, lat, lon, depth, mt,
                      centroid_time, source_time_function)
        events[i] = event
    np.save(root_resources + 'gcmt', events)

def _mt_from_tensor(tensor, Mw):
    m_rr = tensor.m_rr * 1e-18
    m_rt = tensor.m_rt * 1e-18
    m_rp = tensor.m_rp * 1e-18
    m_tt = tensor.m_tt * 1e-18
    m_tp = tensor.m_tp * 1e-18
    m_pp = tensor.m_pp * 1e-18
    # unit conversion. DSM in units of 10**25 [dyne cm]
    mt = MomentTensor(m_rr, m_rt, m_rp, m_tt, m_tp, m_pp, Mw)
    return mt

def read_catalog():
    """Get the GCMT catalog.
    Returns:
        cat (ndarray): ndarray of pydsm.Event objects
    """
    try:
        cat = np.load(root_resources + 'gcmt.npy', allow_pickle=True)
    except:
        print('Dowloading gcmt catalog.\n'
              + 'Takes a few minutes. Done only once.')
        cat = _download_gcmt_catalog()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            convert_catalog(cat)
        cat = np.load(root_resources + 'gcmt.npy', allow_pickle=True)
    return cat

def _download_gcmt_catalog():
    cat = read_events('https://www.ldeo.columbia.edu/~gcmt/projects/CMT/'
                      'catalog/jan76_dec17.ndk')
    start_year = 2018
    end_year = date.today().year
    p=re.compile(r'[a-z]+\d\d\.ndk')
    for year in range(start_year, end_year+1):
        dir = ('https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/'
               'NEW_MONTHLY/' + str(year))
        r = requests.get(dir)
        ndk_files = p.findall(r.text)
        for ndk_file in set(ndk_files):
            try:
                cat_tmp = read_events(dir + '/' + ndk_file)
                cat.extend(cat_tmp.events)
            except:
                pass
    #cat.write(root_resources + 'gcmt.xml', format='quakeml')
    return cat


if __name__ == '__main__':
    cat = read_catalog()
