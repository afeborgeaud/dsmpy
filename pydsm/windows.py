from obspy.taup import TauPyModel
import numpy as np
from pydsm.event import Event
from pydsm.station import Station

class Windows:
    """Time or frequency window.
    
    Args:
        event (Event): Event object with earthquake info
        stations (list(Station)): list of Station objects
        modelName (str): name of a TauP model
        phaseNames (list(str)): list of names of a TauP seismic phase
    """
    def __init__(self, event, stations, modelName, phaseNames):
        self.event = event
        self.stations = stations
        self.modelName = modelName
        self.phaseNames = phaseNames
        self.epicentral_distances = self._get_epicentral_distances()
        self.travel_times = self.compute()
    
    def _get_epicentral_distances(self):
        epicentral_distances = np.zeros(
            len(self.stations), dtype=np.float64)
        for i, station in enumerate(self.stations):
            epicentral_distances[i] = self.event.get_epicentral_distance(
                station)
        return epicentral_distances

    def compute(self):
        taup_model = TauPyModel(self.modelName)
        travel_times = np.zeros(len(self.stations), dtype=np.float64)
        for i, distance in enumerate(self.epicentral_distances):
            arrivals = taup_model.get_travel_times(
                self.event.depth, distance, phase_list=self.phaseNames)
            if len(arrivals) > 0:
                t = arrivals[0].time
            else:
                t = -1.
            travel_times[i] = t
        return travel_times
    
    def get_windows(self, t_before, t_after):
        windows = np.apply_along_axis(
            lambda x: np.array([x-t_before, x+t_after], dtype=np.float64),
            axis=0, arr=self.travel_times)
        return windows
    
    def get_gaussian_windows_in_frequency_domain(
            self, nspc, tlen, window_width):
        """Compute fourier-transformed gaussian windows.
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
        gauss_windows = np.ones((len(self.stations), 2*nspc+1),
                                 dtype=np.complex128)
        tau = tlen / nspc
        for i, t in enumerate(self.travel_times):
            if t == -1:
                continue
            else:
                # add max period / 2 to center the gaussian
                t += tau / 2.
                gauss_windows[i] = (coeff
                    * np.exp(-omegas**2 * window_width**2 / 8 + 1j*omegas*t))
        return gauss_windows
        

if __name__ == '__main__':
    from pydsm.utils.cmtcatalog import read_catalog
    catalog = read_catalog()
    event = Event.event_from_catalog(
    catalog, '200707211534A')
    stations = [
    Station(
        name='FCC', network='CN',
        latitude=58.7592, longitude=-94.0884), ]
    model = 'prem'
    phase = 'ScS'
    windows = Windows(event, stations, model, phase)
    time_windows = windows.get_windows(10, 10)
    print(time_windows)