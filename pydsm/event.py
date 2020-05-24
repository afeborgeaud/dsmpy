import numpy as np
from pydsm.station import Station
from obspy.taup.taup_geo import calc_dist
from pydsm._tish import parameters as tish_parameters
import warnings

class Event:
    """Represent an earthquake point-source.

    Args:
        event_id (str): GCMT event name
        latitude (float): centroid geographic latitude 
            [-90, 90] in degree
        longitude (float): centroid longitude 
            [-180, 180] in degree
        depth (float) centroid depth in km
        mt (ndarray(3, 3)): moment tensor
        source_time_function (SourceTimeFunction): SourceTimeFunction
            object
    
    Attributes:
        event_id (str): GCMT event name
        latitude (float) centroid geographic latitude 
            [-90, 90] in degree
        longitude (float): centroid longitude 
            [-180, 180] in degree
        depth (float) centroid depth in km
        mt (ndarray(3, 3)): moment tensor
        source_time_function (SourceTimeFunction): SourceTimeFunction
            object
    """

    def __init__(self, event_id, latitude, longitude, depth, mt,
                 source_time_function):
        self.event_id = event_id
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.mt = mt
        self.source_time_function = source_time_function

    @classmethod
    def event_from_catalog(cls, cat, event_id):
        """Build Event from GCMT catalog
        Args:
            cat (ndarray): event catalog.
                see pydsm.utils.cmtcatalog.read_catalog()
            event_id (str): GCMT event identifier
                (e.g., '201906291959A')
        Returns:
            event (Event): Event object
        """
        event = None
        try:
            event = cat[cat == event_id][0]
        except:
            warnings.warn('Event {} not found'.format(event_id))
        return event
    
    def get_epicentral_distance(self, station):
        return calc_dist(
                self.latitude, self.longitude,
                station.latitude, station.longitude,
                6371., tish_parameters['flattening'])

    def __repr__(self):
        return self.event_id

    def __eq__(self, event_id):
        return self.event_id == event_id


class MomentTensor:
    """Represent a point-source moment tensor."""

    def __init__(self, Mrr, Mrt, Mrp, Mtt, Mtp, Mpp):
        self.Mrr = Mrr
        self.Mrt = Mrt
        self.Mrp = Mrp
        self.Mtt = Mtt
        self.Mtp = Mtp
        self.Mpp = Mpp
        if np.abs(np.array([Mrr, Mrt, Mrp, Mtt, Mtp, Mpp])).max() > 1e4:
            warnings.warn("Moment tensor should be in units of 10**25 dyne cm")

    def to_array(self):
        mt = np.zeros((3, 3), dtype=np.float64)
        mt[0, 0] = self.Mrr
        mt[0, 1] = self.Mrt
        mt[0, 2] = self.Mrp
        mt[1, 1] = self.Mtt
        mt[1, 2] = self.Mtp
        mt[2, 2] = self.Mpp
        return mt