"""Main module of dsmpy. Contains the classes for core dsmpy objects
and the methods for serial and parallel computation
of synthetic seismograms using the Fortran DSM."""

from dsmpy._tish import _tish, _calthetaphi
from dsmpy._tish import parameters as tish_parameters
from dsmpy._tish import _pinput as _pinput_sh
from dsmpy._tipsv import _pinput, _tipsv
from dsmpy.spc import spctime
from dsmpy.event import Event
from dsmpy.station import Station
from dsmpy.seismicmodel import SeismicModel
from obspy import Trace
from obspy.core.trace import Stats
from obspy.core.util.attribdict import AttribDict
import obspy.signal.filter
import os
import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np
from mpi4py import MPI
import time
import functools
import warnings

default_params = dict(
    re=0.01, ratc=1e-10, ratl=1e-5, omegai=0.0014053864092981234)


def _is_iterable(obj):
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True


class PyDSMInputFile:
    """Input file for dsmpy.

    Args:
        input_file (str): path of dsmpy input file

    """

    def __init__(self, input_file):
        self.input_file = input_file

    def read(self):
        params = dict()
        params['verbose'] = 0
        with open(self.input_file, 'r') as f:
            for line in f:
                if line.strip().startswith('#'):
                    continue
                key, value = self._parse_line(line)
                if key is not None:
                    params[key] = value
        return params

    def _parse_line(self, line):
        key, value = line.strip().split()[:2]
        if key == 'sac_files':
            full_path = os.path.expanduser(value.strip())
            value_parsed = list(glob.iglob(full_path))
        elif key == 'tlen':
            value_parsed = float(value)
        elif key == 'nspc':
            value_parsed = int(value)
        elif key == 'sampling_hz':
            value_parsed = int(value)
        elif key == 'seismic_model':
            value_parsed = value.strip().lower()
        elif key == 'mode':
            value_parsed = int(value)
        elif key == 'output_folder':
            full_path = os.path.expanduser(value.strip())
            value_parsed = full_path
        elif key == 'verbose':
            value_parsed = int(value)
        else:
            print('Warning: key {} undefined. Ignoring.'.format(key))
            return None, None
        return key, value_parsed


class PyDSMOutput:
    """Output from dsmpy compute methods.

    Args:
        spcs (ndarray): seismic spectra computed with DSM.
            Shape is (3, nr, nspc+1).
        stations (ndarray of Station): stations. Shape is (nr,)
        event (Event): seismic event (the source).
        sampling_hz (int): sampling frequency for time-domain waveforms
        tlen (float): duration of the synthetics (in seconds)
            (better to be 2**n/10)
        nspc (int): number of frequency points in the synthetics
            (better to be 2**n)
        omegai (float): numerical damping using in DSM.
            Default is default_params['omegai'].
            Better to leave it at its default value.
        model_id (str): seismic model identifier (e.g., 'PREM', 'mod1').
            (default is None).

    """
    color_count = 0
    colors = (
        'blue', 'red', 'green', 'orange',
        'purple', 'brown', 'pink', 'cyan',
        'magenta', '')

    def __init__(
            self, spcs, stations, event,
            sampling_hz, tlen, nspc,
            omegai=default_params['omegai'], model_id=None):
        self.spcs = spcs
        self.stations = stations
        self.event = event
        self.sampling_hz = sampling_hz
        self.tlen = tlen
        self.nspc = nspc
        self.omegai = omegai
        self.components = ('Z', 'R', 'T')
        self.dt = 1 / self.sampling_hz
        self.us = None
        self.ts = None
        self.model_id = model_id

    @classmethod
    def output_from_pydsm_input(cls, spcs, pydsm_input):
        return cls(
            spcs, pydsm_input.stations, pydsm_input.event,
            pydsm_input.sampling_hz,
            pydsm_input.tlen, pydsm_input.nspc, pydsm_input.omegai)

    def to_time_domain(self, source_time_function=None):
        """Compute time domain waveforms from spetcra.

        Args:
            source_time_function (SourceTimeFunction):
                source time function object (default is None).

        """
        # if self.us is not None:
        #     return

        if source_time_function is None:
            spct = spctime.SpcTime(self.tlen, self.nspc,
                                   self.sampling_hz, self.omegai,
                                   self.event.source_time_function)
        else:
            spct = spctime.SpcTime(self.tlen, self.nspc,
                                   self.sampling_hz, self.omegai,
                                   source_time_function)
        us = spct.spctime(self.spcs)
        self.us = us
        self.ts = np.linspace(0, self.tlen,
                              spct.npts, endpoint=False)

    def set_source_time_function(self, source_time_function):
        """Set the source time function for convolution.

        Args:
            source_time_function (SourceTimeFunction):
                source time function object.

        """
        self.event.source_time_function = source_time_function

    def write(self, root_path, format):
        """Write to file using obspy.io.write.

        Args:
            root_path (str): path of root folder in which to write.
            format (str): output files format ('sac').

        """
        for tr in self.get_traces():
            filename = '.'.join((
                tr.stats.station, tr.stats.network, tr.stats.sac.kevnm,
                tr.stats.component, format))
            tr.write(os.path.join(root_path, filename), format=format)

    def save(self, path):
        """Save self using pickle.dump().

        Args:
            path (str): name of the output file.

        """
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        """Read file into self using pickle.load().

        Args:
            path (str): name of the file that contains self.

        Returns:
            PyDSMOutput: the loaded PyDSMOutput object.

        """
        with open(path, 'rb') as f:
            output = pickle.load(f)
        return output

    def filter(self, freq, freq2=0., type='lowpass', zerophase=False):
        """Filter time-domain waveforms using obspy.signal.filter.

        Args:
            freq (float): filter frequency
            freq2 (float): filter maximum frequency
                (for bandpass filters only)
            type (str): type of filter ('lowpass', 'bandpass').
            zerophase (bool): use zero-phase filter (default is False).

        """
        if self.us is None:
            self.to_time_domain(self.event.source_time_function)

        if type not in {'bandpass', 'lowpass'}:
            raise ValueError('Expect "bandpass" or "lowpass" for arg "type"')

        if type == 'bandpass':
            assert freq2 > freq

        for icomp in range(3):
            for ir in range(self.us.shape[1]):
                if type == 'lowpass':
                    self.us[icomp, ir] = obspy.signal.filter.lowpass(
                        self.us[icomp, ir], freq,
                        df=self.sampling_hz, zerophase=zerophase)
                elif type == 'bandpass':
                    self.us[icomp, ir] = obspy.signal.filter.bandpass(
                        self.us[icomp, ir], freq, freq2,
                        df=self.sampling_hz, zerophase=zerophase)

    def get_traces(self):
        """Return list of obspy.Trace.

        Returns:
            list ot obspy.Trace: traces

        """
        traces = []
        if self.us is None:
            self.to_time_domain()
        for icomp in range(3):
            for ista in range(self.get_nr()):
                station = self.stations[ista]
                data = self.us[icomp, ista]
                stats = Stats()
                stats.network = station.network
                stats.station = station.name
                stats.sampling_rate = self.sampling_hz
                stats.delta = self.dt
                stats.starttime = 0.
                # stats.endtime = self.tlen
                stats.npts = len(data)
                stats.component = self.components[icomp]
                sac_header = AttribDict(**dict(
                    b=0, delta=self.dt, depmax=data.max(),
                    depmin=data.min(), depmen=data.mean(),
                    e=self.tlen, npts=len(data), evdp=self.event.depth,
                    evla=self.event.latitude, evlo=self.event.longitude,
                    kevnm=self.event.event_id, knetwk=station.network,
                    kstnm=station.name, gcarc=0.))
                stats.sac = sac_header
                trace = Trace(data=data, header=stats)
                traces.append(trace)
        return traces

    def window_spcs(self, windows, window_width):
        """Window the spectra in the frequency domain
        using gaussian windows to isolate portions
        of time-domain waveforms.

        Args:
            windows (list of Window): windows.
            window_width (float): duration (in seconds).

        Returns:
            PyDSMOutput: windowed PyDSMOutput object.

        """
        gaussian_windows = [
            window.get_gaussian_window_in_frequency_domain(
                self.nspc, self.tlen, window_width)
            for window in windows]
        spcs_windowed = np.zeros_like(self.spcs)
        for i in range(self.nspc + 1):
            start = self.nspc - i
            end = start + self.nspc + 1
            spcs_windowed[:, :, i] = np.sum(
                self.spcs * gaussian_windows[:, start:end], axis=2)
        output = self.__copy__()
        output.spcs = spcs_windowed
        return output

    def get_nr(self):
        """Returns the number of receivers (stations)."""
        return len(self.stations)

    def _normalize(self, ys, mode='self'):
        if mode == 'self':
            maxs = ys.max(axis=2).reshape(3, -1)
            maxs = np.where(maxs == 0, np.inf, maxs)
            maxs = maxs.reshape((*ys.shape[:2], 1))
            return 0.5 * ys / maxs
        elif mode == 'none':
            # TODO ensure minimum distance
            maxs = ys[:, 0, :].max(axis=1).reshape(3, 1, 1)
            return .8 * ys / maxs

    def _plot(
            self, xs, ys, color='black', axes=None, distance_min=0.,
            distance_max=np.inf, label=None, normalize='self',
            xlabel='Time (s)', slowness=0.):
        if axes is None:
            fig, axes = plt.subplots(3, 1, sharey=True, sharex=True)
        else:
            assert len(axes) == 3
            fig = None
            PyDSMOutput.color_count += 1
        distances = np.zeros(len(self.stations), dtype=np.float64)
        for ir in range(len(self.stations)):
            distances[ir] = self.event.get_epicentral_distance(
                self.stations[ir])
        if distance_min == 0:
            distance_min = distances.min()
        indexes = (distances >= distance_min) & (distances <= distance_max)
        ys_ = self._normalize(ys[:, indexes, :], mode=normalize)
        distances_ = distances[indexes]
        for ir in range(indexes.sum()):
            label_ = label if ir == len(self.stations) - 1 else None
            # reduced time plots
            reduce_time = (distances_[ir] - distance_min) * slowness
            reduce_index = int(reduce_time * self.sampling_hz)
            for icomp in range(3):
                axes[icomp].plot(
                    xs[reduce_index:] - reduce_time,
                    (ys_[icomp, ir, reduce_index:]
                     + distances_[ir]),
                    color=color, label=label_)
                axes[icomp].set_xlabel(xlabel)
                axes[icomp].set_title(self.components[icomp])
                if label is not None:
                    axes[icomp].legend()
            axes[0].set_ylabel('Distance (deg)')
        return fig, axes

    def plot_spc(
            self, color='black', axes=None, distance_min=0.,
            distance_max=np.inf, label=None, normalize='self'):
        """Plot a frequency-domain (spectra) record section.

        Args:
            axes (matplotlib.axes): ax on which to plot
            distance_min (float): minimum epicentral distance (deg)
            distance_max (float): maximum epicentral distance (deg)
            label (str): label for model name
            normalize (str): 'self' for self-normalization
            or 'none' to see amplitude decay with distance

        Returns:
            figure: matplotlib figure object.
            Axes: matplotlib Axes object.

        """
        freqs = np.linspace(
            0, self.nspc / self.tlen, self.nspc + 1, endpoint=True)
        return self._plot(
            freqs, np.abs(self.spcs), color=color, axes=axes,
            distance_min=distance_min,
            distance_max=distance_max, label=label, normalize=normalize,
            xlabel='Frequency (Hz)')

    def plot(
            self, color='black', axes=None, distance_min=0.,
            distance_max=np.inf, label=None, normalize='self',
            slowness=0.):
        if self.us is None:
            self.to_time_domain()
        """Plot a (time-domain) record section.
        
        Args:
            axes (matplotlib.axes): ax on which to plot
            distance_min (float): minimum epicentral distance (deg)
            distance_max (float): maximum epicentral distance (deg)
            label (str): label for model name
            normalize (str): 'self' for self-normalization
                or 'none' to see amplitude decay with distance.
        
        Returns:
            figure: matplotlib figure object.
            Axes: matplotlib Axes object.

        """
        if slowness != 0:
            xlabel = ('Time - {:.1f}*distance (s)'
                      .format(slowness))
        else:
            xlabel = 'Time (s)'
        return self._plot(
            self.ts, self.us, color=color, axes=axes,
            distance_min=distance_min,
            distance_max=distance_max, label=label, normalize=normalize,
            xlabel=xlabel, slowness=slowness)

    def plot_component(
            self, component, windows=None, ax=None,
            align_zero=False, **kwargs):
        """Plot waveforms for one seismic component.

        Args:
            component (Component): seismic component.

        Returns:
            figure: matplotlib figure object.
            Axes: matplotlib Axes object.

        """
        if self.us is None:
            self.to_time_domain()
        if ax is None:
            fig, ax = plt.subplots(1)
        else:
            fig = None
        for i in range(len(self.stations)):
            distance = self.event.get_epicentral_distance(
                self.stations[i])

            # select corresponding window
            if windows is not None:
                windows_tmp = list(filter(
                    lambda w: ((w.station == self.stations[i])
                               and (w.event == self.event)
                               and (w.component == component)),
                    windows))
                if not windows_tmp:
                    continue
                window = windows_tmp[0].to_array()
                if np.isnan(window[0]):
                    data = np.ones(1)
                    ts = np.zeros(1)
                else:
                    i0 = int(window[0] * self.sampling_hz)
                    i1 = int(window[1] * self.sampling_hz)

                    data = self.us[component.value, i, i0:i1]
                    ts = np.linspace(window[0], window[1], len(data))
            else:
                data = self.us[component.value, i]
                ts = np.linspace(
                    0, len(data) / self.sampling_hz, len(data))
            if align_zero:
                ts = np.linspace(
                    0, len(data) / self.sampling_hz, len(data))

            norm = np.abs(data).max() * 2.
            ax.plot(ts, data / norm + distance, **kwargs)
            ax.set(xlabel='Time (s)', ylabel='Distance (deg)')
        return fig, ax

    def __getitem__(self, key):
        """Override __getitem__. Allows following indexations:
        output['Z']
        output['Z', 'sta_net']
        output['Z', ['sta1_net1', 'sta2_net2']]
        """
        if self.us is None:
            self.to_time_domain()
        if sum([type(k) == int for k in key]) == len(key):
            return self.us[key]
        if len(key) == 1:
            if key == 'Z':
                return self.us[0, ...]
            elif key == 'R':
                return self.us[1, ...]
            elif key == 'T':
                return self.us[2, ...]
        elif len(key) == 2:
            if type(key[1]) is int:
                return self.__getitem__(key[0])[key[1]]
            else:
                try:
                    if (type(key[1]) != str) and _is_iterable(key[1]):
                        indexes = [self.stations.index(k)
                                   for k in key[1]]
                        return self.__getitem__(key[0])[indexes, :]
                    else:
                        index = self.stations.index(key[1])
                        return self.__getitem__(key[0])[index, :]
                except:
                    raise KeyError('Station {} not in list'.format(key[1]))
        else:
            raise KeyError('key {} undefined'.format(key))

    def __copy__(self):
        output = PyDSMOutput(
            self.spcs, self.stations, self.event,
            self.sampling_hz, self.tlen, self.nspc, self.omegai)
        return output


class DSMInput:
    """Input parameters used by the Fortran DSM.

    Args:
        re (float):
        ratc (float):
        ratl (float):
        tlen (float): duration of the synthetics (in seconds)
            (must be 2**n/10)
        nspc (int): number of frequency points in the synthetics
            (must be 2**n)
        omegai (float):
        imin (int): 
        imax (int): 
        nzone (int): 
        vrmin (:obj:`ndarray`): 
        vrmax (:obj:`ndarray`): 
        rho (:obj:`ndarray`): 
        vpv (:obj:`ndarray`): 
        vph (:obj:`ndarray`): 
        vsv (:obj:`ndarray`): 
        vsh (:obj:`ndarray`): 
        eta (:obj:`ndarray`): 
        qmu (:obj:`ndarray`): 
        qkappa (:obj:`ndarray`): 
        r0 (float): 
        eqlat (float): 
        eqlon (float): 
        mt (:obj:`ndarray`): 
        nr (int): 
        theta (:obj:`ndarray`): 
        phi (:obj:`ndarray`): 
        lat (float):
        lon (float): 
        output (:obj:`ndarray`): 
        mode (int): 

    """

    def __init__(
            self, re, ratc, ratl, tlen, nspc, omegai, imin, imax, nzone,
            vrmin, vrmax, rho, vpv, vph, vsv, vsh, eta, qmu, qkappa,
            r0, eqlat, eqlon, mt, nr, theta, phi, lat, lon, output, mode=0):
        (self.re, self.ratc,
         self.ratl, self.omegai) = re, ratc, ratl, omegai

        (self.tlen, self.nspc,
         self.imin, self.imax) = tlen, nspc, imin, imax

        (self.nzone, self.vrmin, self.vrmax,
         self.rho, self.vpv, self.vph, self.vsv, self.vsh, self.eta,
         self.qmu, self.qkappa) = (
            nzone, vrmin, vrmax, rho, vpv,
            vph, vsv, vsh, eta, qmu, qkappa)

        self.r0, self.eqlat, self.eqlon, self.mt = (
            r0, eqlat, eqlon, mt)

        (self.nr, self.theta, self.phi, self.lat,
         self.lon, self.output) = (nr, theta, phi,
                                   lat, lon, output)

        # TODO this apparently fixed a bug with zero amp >= 90 degree
        # but should check with Fortran DSM
        self.mt[1, 0] = self.mt[0, 1]
        self.mt[2, 0] = self.mt[0, 2]
        self.mt[2, 1] = self.mt[1, 2]

        self.mode = mode

    def _get_scalar_dict(self):
        return dict(
            re=self.re, ratc=self.ratc, ratl=self.ratl, tlen=self.tlen,
            nspc=self.nspc, omegai=self.omegai, imin=self.imin,
            imax=self.imax, nzone=self.nzone, r0=self.r0,
            eqlat=self.eqlat, eqlon=self.eqlon, nr=self.nr,
            max_nzone=len(self.vrmin), max_nr=len(self.lat))

    @classmethod
    def input_from_file(cls, parameter_file, mode=1):
        """Build a DSMInput object from a DSM input file.

        Args:
            parameter_file (str): path of a DSM input file
            mode (int): 1: P-SV, 2: SH
        Returns:
            DSMInput

        """
        if mode not in {1, 2}:
            raise RuntimeError('mode should be 1 or 2')

        if mode == 1:
            inputs = _pinput(parameter_file)
            (re, ratc, ratl,
             tlen, nspc, omegai,
             imin, imax, nzone) = inputs[:9]
            (vrmin, vrmax, rho,
             vpv, vph, vsv, vsh,
             eta, qmu, qkappa) = inputs[9:19]
            r0, eqlat, eqlon, mt = inputs[19:23]
            nr = inputs[23]
            (theta, phi, lat,
             lon, output) = inputs[24:]
        else:
            inputs = _pinput_sh(parameter_file)
            (re, ratc, ratl,
             tlen, nspc, omegai,
             imin, imax, nzone) = inputs[:9]
            (vrmin, vrmax, rho,
             vsv, vsh, qmu) = inputs[9:15]
            r0, eqlat, eqlon, mt = inputs[15:19]
            nr = inputs[19]
            (theta, phi, lat,
             lon, output) = inputs[20:]
            vpv = None
            vph = None
            eta = None
            qkappa = None

        return cls(
            re, ratc, ratl, tlen, nspc,
            omegai, imin, imax, nzone, vrmin, vrmax,
            rho, vpv, vph, vsv, vsh, eta, qmu, qkappa, r0, eqlat, eqlon,
            mt, nr, theta, phi, lat, lon, output, mode)

    @classmethod
    def input_from_arrays(cls, event, stations,
                          seismic_model, tlen, nspc):
        """Build a DSMInput from user-friendly arguments.
        
        Args:
            event (Event): earthquake information.
            stations (list of Station): seismic stations.
            seismic_model (SeismicModel): Earth structure model 
                (e.g. PREM).
            tlen (float): duration of the synthetics (in seconds)
                (better to be 2**n/10)
            nspc (int): number of frequency points in the synthetics
                (better to be 2**n)

        Returns:
            DSMInput: DSMInput object

        """
        # structure parameters
        # TODO change maxnzone requirements in DSM
        maxnzone = tish_parameters['maxnzone']
        nzone = seismic_model._nzone
        vrmin = np.pad(seismic_model._vrmin, (0, maxnzone - nzone),
                       mode='constant', constant_values=0)
        vrmax = np.pad(seismic_model._vrmax, (0, maxnzone - nzone),
                       mode='constant', constant_values=0)
        rho = np.pad(seismic_model._rho, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        vpv = np.pad(seismic_model._vpv, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        vph = np.pad(seismic_model._vph, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        vsv = np.pad(seismic_model._vsv, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        vsh = np.pad(seismic_model._vsh, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        eta = np.pad(seismic_model._eta, ((0, 0), (0, maxnzone - nzone)),
                     mode='constant', constant_values=0)
        qmu = np.pad(seismic_model._qmu, (0, maxnzone - nzone),
                     mode='constant', constant_values=0)
        qkappa = np.pad(seismic_model._qkappa, (0, maxnzone - nzone),
                        mode='constant', constant_values=0)

        # source parameters
        eqlat = event.latitude
        eqlon = event.longitude
        r0 = 6371. - event.depth
        mt = event.mt.to_array()

        # receiver parameters
        nr = len(stations)
        maxnr = tish_parameters['maxnr']
        lat = np.array([station.latitude for station in stations],
                       dtype=np.float64)
        lon = np.array([station.longitude for station in stations],
                       dtype=np.float64)
        f = functools.partial(_calthetaphi, eqlat=eqlat, eqlon=eqlon)
        theta_phi = [f(stalat, stalon)
                     for stalat, stalon in zip(lat, lon)]
        theta = np.array([x[0] for x in theta_phi])
        phi = np.array([x[1] for x in theta_phi])

        # TODO change maxnr requirements in DSM
        lat = np.pad(lat, (0, maxnr - nr),
                     mode='constant', constant_values=0)
        lon = np.pad(lon, (0, maxnr - nr),
                     mode='constant', constant_values=0)
        theta = np.pad(theta, (0, maxnr - nr),
                       mode='constant', constant_values=0)
        phi = np.pad(phi, (0, maxnr - nr),
                     mode='constant', constant_values=0)

        output = np.empty((tish_parameters['maxnr'], 80), dtype='S1')
        for i, station in enumerate(stations):
            string = (station.name + '_' + station.network
                      + '.' + event.event_id + 'SH.spc')
            arr = np.array([e for e in string], dtype='S1')
            arr = np.pad(
                arr, (0, 80 - len(arr)),
                mode='constant', constant_values='')
            output[i, :] = arr

        # parameters for DSM computation (advanced)
        re = default_params['re']
        ratc = default_params['ratc']
        ratl = default_params['ratl']
        omegai = default_params['omegai']
        imin = 0
        imax = nspc

        return cls(re, ratc, ratl, tlen, nspc, omegai, imin, imax,
                   nzone, vrmin, vrmax, rho, vpv, vph, vsv, vsh,
                   eta, qmu, qkappa, r0, eqlat, eqlon, mt, nr, theta,
                   phi, lat, lon, output, mode=0)

    @classmethod
    def input_from_dict_and_arrays(
            cls, scalar_dict, vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, mt, lat, lon, phi, theta):
        output = np.empty((len(theta), 80), dtype='|S1')
        return cls(
            scalar_dict['re'], scalar_dict['ratc'],
            scalar_dict['ratl'], scalar_dict['tlen'],
            scalar_dict['nspc'], scalar_dict['omegai'],
            scalar_dict['imin'], scalar_dict['imax'],
            scalar_dict['nzone'], vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, scalar_dict['r0'],
            scalar_dict['eqlat'], scalar_dict['eqlon'], mt,
            scalar_dict['nr'], theta, phi, lat, lon, output)

    def get_inputs_for_tish(self):
        # TODO modify fortran? Else, have to take care of case
        # number of core layers != 2
        if self.mode == 0 or self.mode == 1:
            nzone = self.nzone - 2
            vrmin = np.pad(self.vrmin[2:], (0, 2), constant_values=0)
            vrmax = np.pad(self.vrmax[2:], (0, 2), constant_values=0)
            qmu = np.pad(self.qmu[2:], (0, 2), constant_values=0)
            npad = ((0, 0), (0, 2))
            rho = np.pad(self.rho[:, 2:], npad, constant_values=0)
            vsv = np.pad(self.vsv[:, 2:], npad, constant_values=0)
            vsh = np.pad(self.vsh[:, 2:], npad, constant_values=0)
        else:
            nzone = self.nzone
            vrmin = self.vrmin
            vrmax = self.vrmax
            qmu = self.qmu
            rho = self.rho
            vsv = self.vsv
            vsh = self.vsh

        inputs = (self.re, self.ratc, self.ratl, self.tlen,
                  self.nspc, self.omegai, self.imin, self.imax,
                  nzone, vrmin, vrmax, rho,
                  vsv, vsh, qmu, self.r0,
                  self.eqlat, self.eqlon, self.mt, self.nr,
                  self.theta, self.phi, self.lat, self.lon, self.output)
        return inputs

    def get_inputs_for_tipsv(self):
        inputs = (
            self.re, self.ratc, self.ratl, self.tlen,
            self.nspc, self.omegai, self.imin, self.imax,
            self.nzone, self.vrmin, self.vrmax, self.rho,
            self.vpv, self.vph, self.vsv, self.vsh, self.eta,
            self.qmu, self.qkappa, self.r0, self.eqlat, self.eqlon,
            self.mt, self.nr, self.theta, self.phi, self.lat, self.lon,
            self.output)
        return inputs


class PyDSMInput(DSMInput):
    """Input parameters for dsmpy compute methods.

    Args:
        dsm_input (DSMInput): input parameters for Fortran DSM.
        sampling_hz (int): sampling frequency for time-domain waveforms.
        mode (int): 1: P-SV, 2: SH. (default is 1).

    """

    def __init__(
            self, dsm_input, sampling_hz=None,
            mode=1):
        super().__init__(
            *dsm_input.get_inputs_for_tipsv())
        self.sampling_hz = self.find_optimal_sampling_hz(sampling_hz)
        self.stations = self._parse_stations()
        self.event = self._parse_event()
        self.mode = mode

        assert mode in {1, 2}

    @classmethod
    def input_from_file(cls, parameter_file,
                        sampling_hz=None, source_time_function=None,
                        mode=1):
        """Build a PyDSMInput object from a DSM input file.
        
        Args:
            parameter_file (str): path of a DSM input file.
            sampling_hz (float): sampling frequency
            for time-domain waveforms.
            source_time_function (SourceTimeFunction): 
            SourceTimeFunction object.
            mode (int): 1: P-SV, 2: SH.

        Returns:
            PyDSMInput object

        """
        dsm_input = DSMInput.input_from_file(parameter_file, mode)
        pydsm_input = cls(dsm_input, sampling_hz, mode)
        pydsm_input.set_source_time_function(source_time_function)
        return pydsm_input

    @classmethod
    def input_from_arrays(
            cls, event, stations,
            seismic_model, tlen, nspc,
            sampling_hz):
        """Build a PyDSMInput from user-friendly arguments.
        
        Args:
            event (Event): earthquake information.
            stations (list of Station): seismic stations.
            seismic_model (SeismicModel): Earth structure model 
                (e.g. PREM).
            tlen (float): duration of the synthetics (in seconds)
                (better to be 2**n/10)
            nspc (int): number of frequency points in the synthetics
                (better to be 2**n)
            sampling_hz (float): sampling frequency
                for time-domain synthetics.

        Returns:
            PyDSMInput: PyDSMInput object.

        """
        dsm_input = DSMInput.input_from_arrays(event, stations,
                                               seismic_model, tlen, nspc)
        pydsm_input = cls(dsm_input, sampling_hz, mode=1)
        pydsm_input.set_source_time_function(event.source_time_function)
        return pydsm_input

    def find_optimal_sampling_hz(self, sampling_hz):
        if sampling_hz is not None:
            return sampling_hz
        else:
            Tmin = (self.tlen
                    / self.imax)
            optimal_sampling_hz = 40 / Tmin  # TODO: check this scaling
            return optimal_sampling_hz

    def set_source_time_function(self, source_time_function):
        self.event.source_time_function = source_time_function

    def _parse_stations(self):
        stations = []
        for i in range(self.nr):
            name, net = self.output[i].tobytes(). \
                decode('utf-8').split('/')[-1].split('.')[0].split('_')
            station = Station(name, net, self.lat[i], self.lon[i])
            stations.append(station)
        return tuple(stations)

    def _parse_event(self):
        event_id = self.output[0].tobytes(). \
            decode('utf-8').split('/')[-1].split('.')[1]
        if event_id[-2:] == 'SH':
            event_id = event_id[:-2]
        elif event_id[-3:] == 'PSV':
            event_id = event_id[:-3]
        else:
            raise RuntimeError('{}'.format(event_id))
        event = Event(event_id, self.eqlat, self.eqlon,
                      6371. - self.r0, self.mt, None, None)
        return event


def compute(pydsm_input, write_to_file=False,
            mode=0):
    """Compute spectra using DSM.

    Args:
        dsm_input (PyDSMInput): inputs for DSM.
        mode (int): computation mode. 0: both, 1: P-SV, 2: SH.
        write_to_file (bool): write spetrum files to disk
        as specified in dsm_input.output.

    Returns:
        dsm_output (PyDSMOutput): object containing spectra and
        stations/source information.
    
    Note:
        SH and P-SV spectra are summed by default. Using only P-SV
        or SH results in non-physical waves and should be avoided.
        See Kawai et al. (2006) for details.

    """
    if mode not in {0, 1, 2}:
        raise RuntimeError('mode={} undefined. Should be 0, 1, or 2'
                           .format(mode))
    if mode == 0:
        print('compute PSV')
        spcs = _tipsv(
            *pydsm_input.get_inputs_for_tipsv(),
            write_to_file)
        print('compute SH')
        spcs += _tish(
            *pydsm_input.get_inputs_for_tish(),
            write_to_file)
    elif mode == 1:
        print('compute PSV')
        spcs = _tipsv(
            *pydsm_input.get_inputs_for_tipsv(),
            write_to_file)
    else:
        print('compute SH')
        spcs = _tish(
            *pydsm_input.get_inputs_for_tish(),
            write_to_file)

    dsm_output = PyDSMOutput.output_from_pydsm_input(spcs, pydsm_input)
    return dsm_output


def compute_parallel(
        pydsm_input, comm, mode=0, write_to_file=False,
        verbose=0):
    """Compute spectra using DSM with data parallelization.

    """
    rank = comm.Get_rank()
    n_cores = comm.Get_size()

    if rank == 0:
        scalar_dict = pydsm_input._get_scalar_dict()
        scalar_dict['verbose'] = verbose
        scalar_dict['mode'] = mode
    else:
        scalar_dict = None
    scalar_dict = comm.bcast(scalar_dict, root=0)

    if mode not in {0, 1, 2}:
        raise RuntimeError('mode={} undefined. Should be 0, 1, or 2'
                           .format(mode))

    if rank == 0:
        rho = pydsm_input.rho
        vpv = pydsm_input.vpv
        vph = pydsm_input.vph
        vsv = pydsm_input.vsv
        vsh = pydsm_input.vsh
        eta = pydsm_input.eta
        qmu = pydsm_input.qmu
        qkappa = pydsm_input.qkappa
        vrmin = pydsm_input.vrmin
        vrmax = pydsm_input.vrmax
        mt = pydsm_input.mt
    else:
        rho = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        vpv = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        vph = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        vsv = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        vsh = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        eta = np.empty((4, scalar_dict['max_nzone']),
                       dtype=np.float64, order='F')
        qmu = np.empty(scalar_dict['max_nzone'],
                       dtype=np.float64, order='F')
        qkappa = np.empty(scalar_dict['max_nzone'],
                          dtype=np.float64, order='F')
        vrmin = np.empty(scalar_dict['max_nzone'],
                         dtype=np.float64, order='F')
        vrmax = np.empty(scalar_dict['max_nzone'],
                         dtype=np.float64, order='F')
        mt = np.empty((3, 3), dtype=np.float64, order='F')

    comm.Bcast(rho, root=0)
    comm.Bcast(vpv, root=0)
    comm.Bcast(vph, root=0)
    comm.Bcast(vsv, root=0)
    comm.Bcast(vsh, root=0)
    comm.Bcast(eta, root=0)
    comm.Bcast(qmu, root=0)
    comm.Bcast(qkappa, root=0)
    comm.Bcast(vrmin, root=0)
    comm.Bcast(vrmax, root=0)
    comm.Bcast(mt, root=0)

    if rank == 0:
        start_indices, chunk_sizes = _get_chunk_start_indices(
            pydsm_input.nr, n_cores)
        chunk_size = chunk_sizes[0]
        lon = pydsm_input.lon[:pydsm_input.nr]
        lat = pydsm_input.lat[:pydsm_input.nr]
        phi = pydsm_input.phi[:pydsm_input.nr]
        theta = pydsm_input.theta[:pydsm_input.nr]
    else:
        chunk_size = None
        lon = None
        lat = None
        phi = None
        theta = None
        start_indices = None
        chunk_sizes = None
    chunk_size = comm.bcast(chunk_size, root=0)
    if rank == 0:
        comm.send(chunk_sizes[-1], dest=n_cores - 1, tag=11)
    if rank == n_cores - 1:
        chunk_size = comm.recv(source=0, tag=11)

    lon_local = np.empty(chunk_size, dtype=np.float64)
    lat_local = np.empty(chunk_size, dtype=np.float64)
    phi_local = np.empty(chunk_size, dtype=np.float64)
    theta_local = np.empty(chunk_size, dtype=np.float64)

    comm.Scatterv([lon, chunk_sizes, start_indices, MPI.DOUBLE],
                  lon_local, root=0)
    comm.Scatterv([lat, chunk_sizes, start_indices, MPI.DOUBLE],
                  lat_local, root=0)
    comm.Scatterv([phi, chunk_sizes, start_indices, MPI.DOUBLE],
                  phi_local, root=0)
    comm.Scatterv([theta, chunk_sizes, start_indices, MPI.DOUBLE],
                  theta_local, root=0)

    lat_local = np.pad(lat_local, (0, scalar_dict['max_nr'] - chunk_size),
                       mode='constant', constant_values=0)
    lon_local = np.pad(lon_local, (0, scalar_dict['max_nr'] - chunk_size),
                       mode='constant', constant_values=0)
    phi_local = np.pad(phi_local, (0, scalar_dict['max_nr'] - chunk_size),
                       mode='constant', constant_values=0)
    theta_local = np.pad(theta_local, (0, scalar_dict['max_nr'] - chunk_size),
                         mode='constant', constant_values=0)

    scalar_dict['nr'] = chunk_size
    input_local = DSMInput.input_from_dict_and_arrays(
        scalar_dict, vrmin, vrmax, rho, vpv, vph,
        vsv, vsh, eta, qmu, qkappa, mt, lat_local,
        lon_local, phi_local, theta_local)

    start_time = time.time()
    if scalar_dict['mode'] == 0:
        spcs_local = _tish(*input_local.get_inputs_for_tish(),
                           write_to_file=False)
        spcs_local += _tipsv(*input_local.get_inputs_for_tipsv(),
                             write_to_file=False)
    elif scalar_dict['mode'] == 1:
        spcs_local = _tipsv(*input_local.get_inputs_for_tipsv(),
                            write_to_file=False)
    else:
        spcs_local = _tish(*input_local.get_inputs_for_tish(),
                           write_to_file=False)
    end_time = time.time()
    print('{} paths: processor {} in {} s'
          .format(input_local.nr, rank, end_time - start_time))

    # TODO change the order of outputu in DSM 
    # to have nr as the last dimension
    spcs_local = np.array(spcs_local.transpose(0, 2, 1), order='F')

    if rank == 0:
        spcs_gathered = np.empty((3, (pydsm_input.nspc + 1), pydsm_input.nr),
                                 dtype=np.complex128, order='F')
    else:
        spcs_gathered = None

    if rank == 0:
        spcs_chunk_sizes = tuple([3 * size * (scalar_dict['nspc'] + 1)
                                  for size in chunk_sizes])
        spcs_start_indices = tuple([3 * i * (scalar_dict['nspc'] + 1)
                                    for i in start_indices])
    else:
        spcs_chunk_sizes = None
        spcs_start_indices = None

    comm.Barrier()
    comm.Gatherv(spcs_local,
                 [spcs_gathered, spcs_chunk_sizes, spcs_start_indices,
                  MPI.DOUBLE_COMPLEX],
                 root=0)
    if rank == 0:
        spcs_gathered = spcs_gathered.transpose(0, 2, 1)
        output = PyDSMOutput.output_from_pydsm_input(
            spcs_gathered, pydsm_input)
    else:
        output = None

    return output


# TODO implements mode when tipsv ready
# TODO check write_to_file
def compute_dataset_parallel(
        dataset, seismic_model,
        tlen, nspc, sampling_hz,
        comm, mode=0, write_to_file=False,
        verbose=0, log=None):
    """Compute spectra using DSM with data parallelization.

    Args:
        dataset (Dataset): dataset of events & stations.
        comm (MPI.COMM_WORLD): MPI communicator.
        mode (int): computation mode. 0: both, 1: P-SV, 2: SH.
        write_to_file (bool): write output in Kibrary format
        (default: False).
    
    Returns:
        list of PyDSMOutput: list of PyDSMOutput objects with one
            entry for each event in dataset.

    """
    rank = comm.Get_rank()
    n_cores = comm.Get_size()

    if rank == 0:
        scalar_dict = dict(default_params)
        scalar_dict.update(tish_parameters)
        scalar_dict['tlen'] = tlen
        scalar_dict['nspc'] = nspc
        scalar_dict['sampling_hz'] = sampling_hz
        scalar_dict['imin'] = 0
        scalar_dict['imax'] = nspc
        scalar_dict['nzone'] = seismic_model._nzone
        scalar_dict['mode'] = mode
        scalar_dict['verbose'] = verbose
    else:
        scalar_dict = None
    scalar_dict = comm.bcast(scalar_dict, root=0)

    if scalar_dict['mode'] not in {0, 1, 2}:
        raise RuntimeError('mode={} undefined. Should be 0, 1, or 2'
                           .format(mode))

    if rank == 0:
        rho = seismic_model.get_rho()
        vpv = seismic_model.get_vpv()
        vph = seismic_model.get_vph()
        vsv = seismic_model.get_vsv()
        vsh = seismic_model.get_vsh()
        eta = seismic_model.get_eta()
        qmu = seismic_model.get_qmu()
        qkappa = seismic_model.get_qkappa()
        vrmin = seismic_model.get_vrmin()
        vrmax = seismic_model.get_vrmax()
    else:
        rho = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        vpv = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        vph = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        vsv = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        vsh = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        eta = np.empty((4, scalar_dict['maxnzone']),
                       dtype=np.float64, order='F')
        qmu = np.empty(scalar_dict['maxnzone'],
                       dtype=np.float64, order='F')
        qkappa = np.empty(scalar_dict['maxnzone'],
                          dtype=np.float64, order='F')
        vrmin = np.empty(scalar_dict['maxnzone'],
                         dtype=np.float64, order='F')
        vrmax = np.empty(scalar_dict['maxnzone'],
                         dtype=np.float64, order='F')

    comm.Bcast(rho, root=0)
    comm.Bcast(vpv, root=0)
    comm.Bcast(vph, root=0)
    comm.Bcast(vsv, root=0)
    comm.Bcast(vsh, root=0)
    comm.Bcast(eta, root=0)
    comm.Bcast(qmu, root=0)
    comm.Bcast(qkappa, root=0)
    comm.Bcast(vrmin, root=0)
    comm.Bcast(vrmax, root=0)

    if rank == 0:
        sendcounts_sta, displacements_sta = dataset.get_chunks_station(
            n_cores, verbose=verbose)
        sendcounts_eq, displacements_eq = dataset.get_chunks_eq(n_cores)
        sendcounts_mt, displacements_mt = dataset.get_chunks_mt(n_cores)
        lons = dataset.lons
        lats = dataset.lats
        phis = dataset.phis
        thetas = dataset.thetas
        eqlats = dataset.eqlats
        eqlons = dataset.eqlons
        r0s = dataset.r0s
        mts = np.array([mt.to_array() for mt in dataset.mts])
    else:
        sendcounts_sta = None
        displacements_sta = None
        sendcounts_eq = None
        displacements_eq = None
        sendcounts_mt = None
        displacements_mt = None
        lons = None
        lats = None
        phis = None
        thetas = None
        eqlats = None
        eqlons = None
        r0s = None
        mts = None

    nr = np.empty(1, dtype=np.int64)
    comm.Scatter(sendcounts_sta, nr, root=0)

    if scalar_dict['verbose'] >= 1:
        if log:
            log.write('rank {}: nr={}\n'.format(rank, nr))
        else:
            print('rank {}: nr={}'.format(rank, nr))

    lon_local = np.empty(nr, dtype=np.float64)
    lat_local = np.empty(nr, dtype=np.float64)
    phi_local = np.empty(nr, dtype=np.float64)
    theta_local = np.empty(nr, dtype=np.float64)
    eqlon = np.empty((), np.float64)
    eqlat = np.empty((), np.float64)
    r0 = np.empty((), np.float64)
    mt = np.empty((3, 3), np.float64)

    # stations
    comm.Scatterv(
        [lons, sendcounts_sta, displacements_sta, MPI.DOUBLE],
        lon_local, root=0)
    comm.Scatterv(
        [lats, sendcounts_sta, displacements_sta, MPI.DOUBLE],
        lat_local, root=0)
    comm.Scatterv(
        [phis, sendcounts_sta, displacements_sta, MPI.DOUBLE],
        phi_local, root=0)
    comm.Scatterv(
        [thetas, sendcounts_sta, displacements_sta, MPI.DOUBLE],
        theta_local, root=0)
    # events
    comm.Scatterv(
        [eqlats, sendcounts_eq, displacements_eq, MPI.DOUBLE],
        eqlat, root=0)
    comm.Scatterv(
        [eqlons, sendcounts_eq, displacements_eq, MPI.DOUBLE],
        eqlon, root=0)
    comm.Scatterv(
        [r0s, sendcounts_eq, displacements_eq, MPI.DOUBLE],
        r0, root=0)
    # mts
    comm.Scatterv(
        [mts, sendcounts_mt, displacements_mt, MPI.DOUBLE],
        mt, root=0)

    lat_local = np.pad(lat_local, (0, scalar_dict['maxnr'] - nr[0]),
                       mode='constant', constant_values=0)
    lon_local = np.pad(lon_local, (0, scalar_dict['maxnr'] - nr[0]),
                       mode='constant', constant_values=0)
    phi_local = np.pad(phi_local, (0, scalar_dict['maxnr'] - nr[0]),
                       mode='constant', constant_values=0)
    theta_local = np.pad(theta_local, (0, scalar_dict['maxnr'] - nr[0]),
                         mode='constant', constant_values=0)

    scalar_dict['nr'] = nr[0]
    scalar_dict['eqlat'] = eqlat
    scalar_dict['eqlon'] = eqlon
    scalar_dict['r0'] = r0
    input_local = DSMInput.input_from_dict_and_arrays(
        scalar_dict, vrmin, vrmax, rho, vpv, vph,
        vsv, vsh, eta, qmu, qkappa, mt, lat_local,
        lon_local, phi_local, theta_local)

    start_time = time.time()
    if scalar_dict['mode'] == 0:
        spcs_local = _tipsv(
            *input_local.get_inputs_for_tipsv(),
            write_to_file=False)
        spcs_local += _tish(
            *input_local.get_inputs_for_tish(),
            write_to_file=False)
    elif scalar_dict['mode'] == 1:
        spcs_local = _tipsv(
            *input_local.get_inputs_for_tipsv(),
            write_to_file=False)
    else:
        spcs_local = _tish(
            *input_local.get_inputs_for_tish(),
            write_to_file=False)
    end_time = time.time()
    if log:
        log.write('rank {}: {} paths finished in {} s\n'
                  .format(rank, input_local.nr, end_time - start_time))
    else:
        print('rank {}: {} paths finished in {} s'
              .format(rank, input_local.nr, end_time - start_time))

    # TODO change the order of outputu in DSM 
    # to have nr as the last dimension
    spcs_local = np.array(spcs_local.transpose(0, 2, 1), order='F')

    if rank == 0:
        nspc = scalar_dict['nspc']
        spcs_gathered = np.empty((3, nspc + 1, dataset.nr),
                                 dtype=np.complex128, order='F')
    else:
        spcs_gathered = None

    if rank == 0:
        counts_spcs = tuple([3 * size * (nspc + 1)
                             for size in sendcounts_sta])
        displacements_spcs = tuple([3 * i * (nspc + 1)
                                    for i in displacements_sta])
    else:
        counts_spcs = None
        displacements_spcs = None

    comm.Barrier()
    comm.Gatherv(spcs_local,
                 [spcs_gathered, counts_spcs, displacements_spcs,
                  MPI.DOUBLE_COMPLEX],
                 root=0)

    if rank == 0:
        splits = dataset.nrs.cumsum()
        splits = np.concatenate([[0], splits])
        outputs = []
        for i in range(len(splits) - 1):
            start, end = splits[i], splits[i + 1]
            output = PyDSMOutput(
                spcs_gathered[:, :, start:end].transpose(0, 2, 1),
                dataset.stations[start:end],
                dataset.events[i],
                scalar_dict['sampling_hz'],
                scalar_dict['tlen'],
                scalar_dict['nspc'],
                scalar_dict['omegai'])
            outputs.append(output)
    else:
        outputs = None
    comm.Barrier()

    return outputs


def _get_models_array(models, maxnzone):
    """
    Args:
        models (:obj:`list` of :obj:`SeismicModel`):

    Returns:

    """
    n = len(models)
    model_arr = np.empty((4, maxnzone, 6, n), dtype=np.float64, order='F')
    model_scal_arr = np.empty((maxnzone, 4, n), dtype=np.float64, order='F')
    model_ids = np.empty((20, n), dtype='U', order='F')

    for i in range(n):
        model_arr[:, :, 0, i] = models[i].get_rho()
        model_arr[:, :, 1, i] = models[i].get_vpv()
        model_arr[:, :, 2, i] = models[i].get_vph()
        model_arr[:, :, 3, i] = models[i].get_vsv()
        model_arr[:, :, 4, i] = models[i].get_vsh()
        model_arr[:, :, 5, i] = models[i].get_eta()
        model_scal_arr[:, 0, i] = models[i].get_qmu()
        model_scal_arr[:, 1, i] = models[i].get_qkappa()
        model_scal_arr[:, 2, i] = models[i].get_vrmin()
        model_scal_arr[:, 3, i] = models[i].get_vrmax()
        mid = models[i]._model_id
        model_ids[:len(mid), i] = [c for c in mid]
    return model_arr, model_scal_arr, model_ids


def compute_models_parallel(
        dataset, models,
        tlen, nspc, sampling_hz,
        comm, mode=0, write_to_file=False,
        verbose=0):
    """Perform a model grid search with model parallelization.

    Args:
        dataset (Dataset): dataset.
        models (list of SeismicModel):
            list of seismic models.
        tlen (float): duration of the synthetics (in seconds)
            (better to be 2**n/10)
        nspc (int): number of frequency points in the synthetics
            (better to be 2**n)
        sampling_hz (int): 
        mode (int): computation mode.
            0: both, 1: P-SV, 2: SH (default is 0).
        write_to_file (bool): write output in Kibrary format
            (default is False).
        verbose (int): debugging parameter (default is 0).

    Returns:
       list of list of PyDSMOutput: Shape is (n_models, n_events).

    """
    rank = comm.Get_rank()
    n_cores = comm.Get_size()

    # TODO bad bug fix for when len(models) % n_cores != 0
    if rank == 0:
        if len(models) % n_cores != 0:
            n_comp = (int(len(models) / n_cores) + 1) * n_cores - len(models)
            models_ = models + models[:n_comp]
            warnings.warn(
                "n_models % n_cores != 0. "
                "Computing padded models with len={}".format(len(models_)))
        else:
            models_ = models
    else:
        models_ = None

    maxnzone = tish_parameters['maxnzone']

    if rank == 0:
        model_indexes = np.array(list(range(len(models_))), dtype='i')
        n_models = len(models_) // n_cores
        n0_models = len(models_) - n_models * (n_cores - 1)
        sendcounts = np.array([n_models for i in range(n_cores)], dtype='i')
        sendcounts[0] = n0_models
        displacements = sendcounts.cumsum() - sendcounts[0]
    else:
        sendcounts = None
        displacements = None
        model_indexes = None

    nmod = np.empty(1, dtype='i')
    comm.Scatter(sendcounts, nmod, root=0)
    nmod = int(nmod[0])

    model_indexes_local = np.empty(nmod, dtype='i')

    if verbose >= 1:
        print('rank {}: nmod={}'.format(rank, nmod))

    comm.Scatterv(
        [model_indexes, sendcounts, displacements, MPI.INT],
        model_indexes_local, root=0)

    # broadcast models
    if rank == 0:
        model_arr, model_scal_arr, model_ids = (
            _get_models_array(models_, maxnzone))
        sendcounts_mod = tuple(
            [size * 4 * maxnzone * 6 for size in sendcounts])
        displacements_mod = tuple(
            [i * 4 * maxnzone * 6 for i in displacements])
        sendcounts_scal_mod = tuple(
            [size * maxnzone * 4 for size in sendcounts])
        displacements_scal_mod = tuple(
            [i * maxnzone * 4 for i in displacements])
        sendcounts_id = tuple([size * 20 for size in sendcounts])
        displacements_id = tuple([i * 20 for i in displacements])
    else:
        model_arr, model_scal_arr = None, None
        sendcounts_mod = None
        displacements_mod = None
        sendcounts_scal_mod = None
        displacements_scal_mod = None
        sendcounts_id = None
        displacements_id = None
        model_ids = None

    model_arr_local = np.empty(
        (4, maxnzone, 6, nmod), dtype=np.float64, order='F')
    model_scal_arr_local = np.empty(
        (maxnzone, 4, nmod), dtype=np.float64, order='F')
    model_ids_local = np.empty((20, nmod), dtype='U', order='F')

    comm.Scatterv(
        [model_arr, sendcounts_mod, displacements_mod, MPI.DOUBLE],
        model_arr_local, root=0)
    comm.Scatterv(
        [model_scal_arr, sendcounts_scal_mod,
         displacements_scal_mod, MPI.DOUBLE],
        model_scal_arr_local, root=0)
    # comm.Scatterv(
    #     [model_ids, sendcounts_id,
    #     displacements_id, MPI.CHAR],
    #     model_ids_local, root=0)

    # TODO broadcast dataset

    model_event_spc_local = np.empty(
        (len(dataset.events), 3, nspc + 1, dataset.nr, nmod),
        dtype=np.complex128, order='F')
    countmod = 0
    for imod in range(nmod):
        nzone = np.where(model_scal_arr_local[:, 3, imod] == 0)[0][0]
        model_i = SeismicModel(
            model_scal_arr_local[:nzone, 2, imod],
            model_scal_arr_local[:nzone, 3, imod],
            model_arr_local[:, :nzone, 0, imod],
            model_arr_local[:, :nzone, 1, imod],
            model_arr_local[:, :nzone, 2, imod],
            model_arr_local[:, :nzone, 3, imod],
            model_arr_local[:, :nzone, 4, imod],
            model_arr_local[:, :nzone, 5, imod],
            model_scal_arr_local[:nzone, 0, imod],
            model_scal_arr_local[:nzone, 1, imod],
            None, None, None)
        for iev in range(len(dataset.events)):
            start, end = dataset.get_bounds_from_event_index(iev)
            input_local = PyDSMInput.input_from_arrays(
                dataset.events[iev], dataset.stations[start:end],
                model_i, tlen, nspc, sampling_hz)

            if mode == 0:
                spcs_local = _tipsv(
                    *input_local.get_inputs_for_tipsv(),
                    write_to_file=False)
                spcs_local += _tish(
                    *input_local.get_inputs_for_tish(),
                    write_to_file=False)
            elif mode == 1:
                spcs_local = _tipsv(
                    *input_local.get_inputs_for_tipsv(),
                    write_to_file=False)
            else:
                spcs_local = _tish(
                    *input_local.get_inputs_for_tish(),
                    write_to_file=False)
            # TODO change the order of outputu in DSM 
            # to have nr as the last dimension
            spcs_local = np.array(spcs_local.transpose(0, 2, 1), order='F')
            model_event_spc_local[iev, :, :, start:end, countmod] = spcs_local
            if verbose == 2:
                print(rank, model_event_spc_local[0, 2, 15, 0, 0])
            if verbose >= 1:
                # TODO delete some content
                if np.isnan(spcs_local.any()):
                    print('{} some spc is NaN'.format(rank))
                if np.any(spcs_local[2, 1:, :] == 0):
                    mask = (spcs_local[2, :, :] == 0).all(axis=0)
                    print('spc', spcs_local[2, :, mask])
                    print('vsh', model_arr_local[:, :nzone, 4, imod])
                    raise Exception('{} {} some spc is 0'.format(rank, imod))
        countmod += 1

    comm.Barrier()
    if rank == 0:
        counts_spcs = tuple(
            [size * len(dataset.events) * 3 * (nspc + 1) * dataset.nr
             for size in sendcounts])
        displacements_spcs = tuple(
            [j * len(dataset.events) * 3 * (nspc + 1) * dataset.nr
             for j in displacements])

        model_event_spc_gather = np.empty(
            (len(dataset.events), 3, nspc + 1, dataset.nr, len(models_)),
            dtype=np.complex128, order='F')
    else:
        model_event_spc_gather = None
        counts_spcs = None
        displacements_spcs = None

    comm.Gatherv(
        model_event_spc_local,
        [model_event_spc_gather, counts_spcs, displacements_spcs,
         MPI.DOUBLE_COMPLEX], root=0)

    if rank == 0:
        model_event_spc_gather = model_event_spc_gather.transpose(4, 0, 1, 3,
                                                                  2)
        outputs = list()
        for imod in range(len(models)):
            output_event_list = list()
            for iev in range(len(dataset.events)):
                # (n_mod, n_events, 3, nspc+1, n_sta
                start, end = dataset.get_bounds_from_event_index(iev)
                output = PyDSMOutput(
                    model_event_spc_gather[imod, iev],
                    dataset.stations[start:end],
                    dataset.events[iev],
                    sampling_hz,
                    tlen,
                    nspc)
                output_event_list.append(output)
            outputs.append(output_event_list)
    else:
        outputs = None

    return outputs


def _get_chunk_start_indices(nr, n_cores):
    chunk_size = nr // n_cores
    start_indices = [i * chunk_size for i in range(n_cores)]
    chunk_sizes = [chunk_size for i in range(n_cores - 1)]
    chunk_sizes += [nr - start_indices[-1]]
    return tuple(start_indices), tuple(chunk_sizes)
