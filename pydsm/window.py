from obspy.taup import TauPyModel
import numpy as np
from pydsm.event import Event
from pydsm.station import Station
from pydsm.component import Component

class Window:
    """Time or frequency window.
    
    Args:
        travel_time (float): arrival time
        event (pydsm.Event): seismic event
        station (pydsm.Station): seismic station
        phaseName (str): a TauP seismic phase
        component (pydsm.Component): seismic component
        t_before (float): time before arrival (default: 10)
        t_after (float): time after arrival (default: 40)
    """
    def __init__(
            self, travel_time, event, station,
            phase_name, component, t_before=10., t_after=40.):
        self.travel_time = travel_time
        self.event = event
        self.station = station
        self.phase_name = phase_name
        self.component = component
        self.t_before = t_before
        self.t_after = t_after

    def get_epicentral_distance(self):
        '''Returns the epicentral distance in degree.'''
        return self.event.get_epicentral_distance(
                self.station)
    
    def to_array(self):
        '''Returns an ndarray [t_start, t_end].'''
        return np.array(
            [self.travel_time-self.t_before, self.travel_time+self.t_after],
             dtype=np.float64)
    
    def get_gaussian_window_in_frequency_domain(
            self, nspc, tlen, window_width):
        """Compute a gaussian window in the frequency domain.
        Args:
            nspc (int): number of points in frequency domain
            tlen (float): duration of synthetics (in seconds)
            window_width (float): gaussian width (in seconds) = 2*sigma
            omega_shift (float): omega shift
        Returns:
            windows (list(ndarray)): list of gaussian windows
        """
        omega_start = -2 * np.pi * nspc / tlen
        omega_end = -omega_start
        omegas = np.linspace(omega_start, omega_end, 2*nspc+1, endpoint=True)
        coeff = np.sqrt(0.5 * np.pi) * window_width
        gauss_window = np.ones(2*nspc+1,
                                 dtype=np.complex128)
        tau = tlen / nspc
        t = self.travel_time
        if t == np.NaN:
            return None
        else:
            # add max period / 2 to center the gaussian
            t += tau / 2.
            gauss_window = (coeff
                * np.exp(-omegas**2 * window_width**2 / 8 + 1j*omegas*t))
        return gauss_window

    def __repr__(self):
        return '{} {} {:.2f} {} {}'.format(
            self.event, self.station, self.travel_time,
            self.phase_name, self.component.name)
        