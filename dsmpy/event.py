from dsmpy.station import Station
from dsmpy._tish import parameters as tish_parameters
from obspy.taup.taup_geo import calc_dist, calc_dist_azi
import warnings
import numpy as np

class Event:
    """Represent an earthquake point-source.

    Args:
        event_id (str): GCMT event name
        latitude (float): centroid geographic latitude 
            [-90, 90] in degree
        longitude (float): centroid longitude 
            [-180, 180] in degree
        depth (float) centroid depth in km
        mt (MomentTensor): moment tensor
        centroid_time (DateTime): centroid time
        source_time_function (SourceTimeFunction): SourceTimeFunction
        object
    
    Attributes:
        event_id (str): GCMT event name
        latitude (float) centroid geographic latitude 
        [-90, 90] in degree
        longitude (float): centroid longitude 
        [-180, 180] in degree
        depth (float) centroid depth in km
        mt (MomentTensor): moment tensor
        source_time_function (SourceTimeFunction): SourceTimeFunction
        object
        centroid_time (datetime): centroid time.

    """

    def __init__(self, event_id, latitude, longitude, depth, mt,
                 centroid_time, source_time_function):
        self.event_id = event_id
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.mt = mt
        self.source_time_function = source_time_function
        self.centroid_time = centroid_time

    @classmethod
    def event_from_catalog(cls, cat, event_id):
        """Build Event from GCMT catalog.

        Args:
            cat (:obj:`ndarray`): event catalog.
            see pydsm.utils.cmtcatalog.read_catalog().
            event_id (str): GCMT event identifier
            (e.g., '201906291959A').

        Returns:
            event (:obj:`Event`): Event object.

        """
        event = None
        try:
            event = cat[cat == event_id][0]
        except:
            warnings.warn('Event {} not found'.format(event_id))
        return event
    
    def get_epicentral_distance(self, station):
        """Returns the epicentral distance in degrees.

        Args:
            station (Station): station

        Returns:
            float: epicentral distance in degrees

        """
        return calc_dist(
                self.latitude, self.longitude,
                station.latitude, station.longitude,
                6371., tish_parameters['flattening'])

    def get_epicentral_distance_(self, sta_lat, sta_lon):
        """Returns the epicentral distance in degrees.

        Args:
            sta_lat (float): station latitude in degrees
            sta_lon (float): station longitude in degrees

        Returns:
            float: epicentral distance in degrees

        """
        return calc_dist(
                self.latitude, self.longitude,
                sta_lat, sta_lon,
                6371., tish_parameters['flattening'])

    def get_azimuth(self, station):
        """Returns the source-station azimuth in degrees.

        Args:
            station (Station): seismic station

        Returns:
            float: azimuth in degrees

        """
        dist, az, backaz = calc_dist_azi(
            self.latitude, self.longitude,
            station.latitude, station.longitude,
            6371., tish_parameters['flattening'])
        return az

    def get_backazimuth(self, station):
        """Returns the station-source backazimuth in degrees.

        Args:
            station (Station): seismic station

        Returns:
            float: backazimuth in degrees

        """
        dist, az, backaz = calc_dist_azi(
            self.latitude, self.longitude,
            station.latitude, station.longitude,
            6371., tish_parameters['flattening'])
        return backaz

    def __repr__(self):
        return self.event_id

    def __eq__(self, event_id):
        return self.event_id == event_id

    def __hash__(self):
        return hash(self.event_id)


class MomentTensor:
    """Represent a point-source moment tensor.
    
    """

    def __init__(self, Mrr, Mrt, Mrp, Mtt, Mtp, Mpp, Mw=None):
        self.Mrr = Mrr
        self.Mrt = Mrt
        self.Mrp = Mrp
        self.Mtt = Mtt
        self.Mtp = Mtp
        self.Mpp = Mpp
        self.Mw = Mw
        if np.abs(np.array([Mrr, Mrt, Mrp, Mtt, Mtp, Mpp])).max() > 1e4:
            warnings.warn("Moment tensor should be in units of 10**25 dyne cm")

    @classmethod
    def from_dsm_array(cls, mt_arr):
        '''Create a MomentTensor from the mt array from DSM pinput
        (assume the DSM order for mt component).

        Args:
            mt_arr (:obj:`ndarray`): moment tensor array

        Returns:
            mt (:obj:`MomentTensor`): moment tensor
            
        '''
        return cls(
            mt_arr[0,0], mt_arr[0,1], mt_arr[0,2],
            mt_arr[1,1], mt_arr[1,2], mt_arr[2,2])

    def to_array(self):
        mt = np.zeros((3, 3), dtype=np.float64)
        mt[0, 0] = self.Mrr
        mt[0, 1] = self.Mrt
        mt[0, 2] = self.Mrp
        mt[1, 1] = self.Mtt
        mt[1, 2] = self.Mtp
        mt[2, 2] = self.Mpp
        
        # TODO this apparently fixed a bug with zero amp >= 90 degree
        # but should check with Fortran DSM
        mt[1, 0] = mt[0, 1]
        mt[2, 0] = mt[0, 2]
        mt[2, 1] = mt[1, 2]
        return mt
