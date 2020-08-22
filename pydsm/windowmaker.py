from obspy.taup import TauPyModel
import numpy as np
from pydsm.event import Event
from pydsm.station import Station
from pydsm.window import Window
from pydsm.component import Component

class WindowMaker:
    """Utility class to compute list of pydsm.Windows.
    """

    @staticmethod
    def windows_from_obspy_traces(
            trace, model_name, phase_names,
            t_before=10., t_after=40.):
        event = Event(
            trace.stats.sac.kevnm,
            trace.stats.sac.evla,
            trace.stats.sac.evlo,
            trace.stats.sac.evdp,
            None, None, None)
        station = Station(
            trace.stats.station,
            trace.stats.network,
            trace.stats.sac.stla,
            trace.stats.sac.stlo)
        component = Component.parse_component(trace.stats.sac.kcmpnm)

        windows = WindowMaker.compute(
            event, [station], model_name, phase_names,
            [component], t_before, t_after)
        return windows

    @staticmethod
    def windows_from_dataset(
            dataset, model_name, phase_names,
            components, t_before=10., t_after=40.):
        '''Compute windows from a pydsm.Dataset.
        Args:
            dataset (pydsm.Dataset): dataset
            model_name (str): name of the 1-D reference model
            phase_names (list(str)): name of seismic phases
            components (list(pydsm.Component)): seismic components
            t_before (float): time before arrival (default: 10)
            t_after (float): time after arrival (default: 40)
        Returns:
            windows (list(pydsm.Windows)): time windows
        '''
        windows = []
        for i, event in enumerate(dataset.events):
            start, end = dataset.get_bounds_from_event_index(i)
            stations = dataset.stations[start:end]
            tmp_windows = WindowMaker.compute(
                event, stations, model_name, phase_names,
                components, t_before, t_after)
            windows += tmp_windows
        return windows

    @staticmethod
    def compute(
            event, stations, model_name, phase_names,
            components, t_before=10., t_after=40.):
        '''Compute time windows using TauP.
        Args:
            event (pydsm.Event): seismic event
            stations (list(pydsm.Station)): seismic stations
            model_name (str): name of reference 1-D model
            phase_names (list(str)): list of TauP phase names
            components (list(pydsm.Component)): seismic components
            t_before (float): time before arrival (default: 10)
            t_after (float): time after arrival (default: 40)
        Returns:
            windows (list(pydsm.Window)): time windows
        '''
        taup_model = TauPyModel(model_name)
        windows = []
        for station in stations:
            distance = event.get_epicentral_distance(station)
            arrivals = taup_model.get_travel_times(
                event.depth, distance, phase_list=phase_names)
            # at the moment consider only first arrival
            processed_phases = set()
            if len(arrivals) > 0:
                for arrival in arrivals:
                    if arrival.name not in processed_phases:
                        for component in components:
                            windows.append(
                                Window(
                                    arrival.time, event, station, arrival.name,
                                    component, t_before, t_after))
                    processed_phases.add(arrival.name)
            else:
                for component in components:
                    windows.append(Window(
                        np.NaN, event, station, None,
                        component, np.NaN, np.NaN))
        return windows

    @staticmethod
    def set_limit(windows, t_before, t_after):
        '''Set t_before and t_after for all window in windows.'''
        for i in range(len(windows)):
            windows[i].t_before = t_before
            windows[i].t_after = t_after
    
    # def get_window_array(windows, t_before, t_after):
    #     '''Return a ndarray of [t_start, t_end].
    #     Args:
    #         windows (list(pydsm.Window)): time windows
    #         t_before (float):
    #         t_after (float):
    #     Returns:
    #         window_array (ndarra):
    #     '''
    #     windows = np.apply_along_axis(
    #         lambda x: np.array([x-t_before, x+t_after], dtype=np.float64),
    #         axis=0, arr=self.travel_times)
    #     return windows.T
    
    # def get_gaussian_windows_in_frequency_domain(
    #         self, nspc, tlen, window_width):
    #     """Compute fourier-transformed gaussian windows.
    #     Args:
    #         nspc (int): number of points in frequency domain
    #         tlen (float): duration of synthetics (in seconds)
    #         window_width (float): gaussian width (in seconds) = 2*sigma
    #         omega_shift (float): omega shift
    #     Returns:
    #         windows (list(ndarray)): list of gaussian windows
    #     """
    #     omega_start = -2 * np.pi * nspc / tlen
    #     omega_end = -omega_start
    #     omegas = np.linspace(omega_start, omega_end, 2*nspc+1, endpoint=True)
    #     coeff = np.sqrt(0.5 * np.pi) * window_width
    #     gauss_windows = np.ones((len(self.stations), 2*nspc+1),
    #                              dtype=np.complex128)
    #     tau = tlen / nspc
    #     for i, t in enumerate(self.travel_times):
    #         if t == np.NaN:
    #             continue
    #         else:
    #             # add max period / 2 to center the gaussian
    #             t += tau / 2.
    #             gauss_windows[i] = (coeff
    #                 * np.exp(-omegas**2 * window_width**2 / 8 + 1j*omegas*t))
    #     return gauss_windows
        

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
    phases = ['S', 'ScS']
    components = [Component.T]
    windows = WindowMaker.compute(event, stations, model, phases, components)
    print(windows)