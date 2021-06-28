from dsmpy.dsm import PyDSMInput
from dsmpy.event import Event, MomentTensor
from dsmpy.station import Station
from dsmpy._tish import _calthetaphi
from dsmpy import root_resources
from dsmpy.utils.cmtcatalog import read_catalog
from dsmpy.component import Component
from dsmpy.spc.stf import SourceTimeFunction
from dsmpy.spc.stfcatalog import STFCatalog
import numpy as np
from obspy import read
import obspy.signal.filter
import pandas as pd
import matplotlib.pyplot as plt
from mpi4py import MPI
from collections import defaultdict

DATA_FLOAT_PREC = np.float32

class Dataset:
    """Represents a dataset of events and stations.

    The data array is not None only if the dataset was defined using
    Dataset.read_from_sac(headonly=False). In this case,
    the data array is of shape (1, 3, n_records, npts),
    where n_records is the number of seismic records,
    or event-station pairs, and npts is the number of time points
    for the longest record.
    Dimension 1 corresponds to the 3 seismic components (Z, R, T).
    Dimension 0 has length >= 1 only after dataset.apply_windows().
    In this case, dimension 0 encodes the number of time windows (i.e.,
    the number of different phases).
    
    Args:
        lats (ndarray): stations latitudes for each record (nr,).
        lons (ndarray): stations longitudes for each record (nr,).
        phis (ndarray): stations phis for each record (nr,).
        thetas (ndarray): stations thetas for each record (nr,).
        eqlats (ndarray): centroids latitudes (nev,).
        eqlons (ndarray): centroids longitudes (nev,).
        r0s (ndarray): centroids radii (nev,).
        mts (ndarray of MomentTensor): array of moment tensors (nev,).
        nrs (ndarray of int): number of stations for each event (nev,).
        nr (int): total number of event-station pairs.
        stations (ndarray of Station): seismic stations (nr,).
        events (ndarray of Event)): seismic events (nev,).
        data (ndarray): 3-components waveform data.
        nw: number of windows used to cut the data (nw,3,nr,npts).
            If self.cut_data() hasn't been called, then nw=1.
        sampling_hz (int): sampling frequency for data.
            Used for computation with pydsm.
        
    """
    def __init__(
            self, lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events, data=None,
            sampling_hz=20, is_cut=False):
        self.lats = lats
        self.lons = lons
        self.phis = phis
        self.thetas = thetas
        self.eqlats = eqlats
        self.eqlons = eqlons
        self.r0s = r0s
        self.mts = mts
        self.nrs = nrs
        self.nr = len(self.lats)
        self.stations = stations
        self.events = events
        if data is not None:
            self.data = np.array(data)
        else:
            self.data = None
        self.sampling_hz = sampling_hz
        self.is_cut = is_cut

    def copy(self):
        """Return a deep copy of self.

        Returns:
            Dataset: deep copy of self

        """
        return Dataset(
            self.lats, self.lons, self.phis, self.thetas,
            self.eqlats, self.eqlons, self.r0s, self.mts,
            self.nrs, self.stations, self.events, self.data,
            self.sampling_hz, self.is_cut)

    @classmethod
    def dataset_from_files(cls, parameter_files, file_mode=1):
        """Create a Dataset object from a list of DSM input files.
        This dataset does not contain waveform data (self.data is None),
        and is used only to compute synthetics.

        Args:
            parameter_file (str): path of a DSM input file.
            file_mode (int): The kind of DSM input file. 1: P-SV, 2: SH.

        Returns:
            Dataset

        """
        pydsm_inputs = [PyDSMInput.input_from_file(file, file_mode=file_mode)
                            for file in parameter_files]
        
        lats = np.concatenate([input.lat[:input.nr]
                               for input in pydsm_inputs])
        lons = np.concatenate([input.lon[:input.nr]
                               for input in pydsm_inputs])
        phis = np.concatenate([input.phi[:input.nr]
                               for input in pydsm_inputs])
        thetas = np.concatenate([input.theta[:input.nr]
                                for input in pydsm_inputs])
        eqlats = np.array([input.eqlat for input in pydsm_inputs])
        eqlons = np.array([input.eqlon for input in pydsm_inputs])
        r0s = np.array([input.r0 for input in pydsm_inputs])
        mts = np.array(
            [MomentTensor.from_dsm_array(input.mt)
            for input in pydsm_inputs])

        nrs = np.array([input.nr for input in pydsm_inputs])

        stations = np.concatenate([input.stations
                                        for input in pydsm_inputs])
        events = np.array([input.event
                                for input in pydsm_inputs])

        return cls(
            lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events)
    
    @classmethod
    def dataset_from_arrays(cls, events, stations, sampling_hz=20):
        """Create a Dataset object from a list of events and stations.
        This dataset does not contain waveform data (self.data is None),
        and is used only to compute synthetics.

        Args:
            events (iterable of Event): earthquake events
            stations (iterable of Station): seismic stations
            sampling_hz (float): waveform sampling that will be
                inherited by the synthetics (default is 20)

        Returns:
            Dataset
        """
        eqlats = np.array([e.latitude for e in events])
        eqlons = np.array([e.longitude for e in events])
        r0s = np.array([6371. - e.depth for e in events])
        mts = np.array([e.mt for e in events])
        nrs = np.array([len(stas) for stas in stations])
        lats = np.array([s.latitude for stas in stations for s in stas])
        lons = np.array([s.longitude for stas in stations for s in stas])

        thetas = np.zeros(len(lats), dtype=np.float64)
        phis = np.zeros(len(lats), dtype=np.float64)
        count = 0
        for i in range(len(events)):
            for sta in stations[i]:
                theta_phi = _calthetaphi(
                    sta.latitude, sta.longitude,
                    events[i].latitude, events[i].longitude)
                thetas[count], phis[count] = theta_phi
                count += 1
        
        stations_ = np.concatenate(stations)
        events_ = np.array(events)

        return cls(lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations_, events_, data=None,
            sampling_hz=sampling_hz)


    @classmethod
    def dataset_from_sac(
            cls, sac_files, verbose=0, headonly=True):
        """Creates a dataset from a list of sac files.
        With headonly=False, time series data from the sac_files
        will be stored in self.data.

        For parallel applications using MPI, headonly=False (i.e.,
        reading the data from sac files) only applies to rank 0, so
        as not to saturate the memory.

        Args:
            sac_files (list of str): list of paths to sac files.
            verbose (int): 0: quiet, 1: debug.
            headonly (bool): if True, read only the metadata.
                If False, includes data.

        Returns:
            Dataset: dataset

        Examples:
            >>> sac_files = ['FCC.CN.201205280204A.T']
            >>> dataset = Dataset.dataset_from_sac(
            ...        sac_files, headonly=False)

        """
        if MPI.COMM_WORLD.Get_rank() == 0:
            traces = [read(sac_file, headonly=headonly)[0]
                      for sac_file in sac_files]
        else:
            traces = [read(sac_file, headonly=True)[0]
                      for sac_file in sac_files]

        sampling_hz = int(traces[0].stats.sampling_rate)

        lats_ = []
        lons_ = []
        names_ = []
        nets_ = []
        eqlats_ = []
        eqlons_ = []
        eqdeps_ = []
        evids_ = []
        data_ = []
        components_ = []
        indices_ = list(range(len(traces)))

        for tr in traces:
            lats_.append(tr.stats.sac.stla)
            lons_.append(tr.stats.sac.stlo)
            names_.append(tr.stats.sac.kstnm)
            nets_.append(tr.stats.sac.knetwk)
            eqlats_.append(tr.stats.sac.evla)
            eqlons_.append(tr.stats.sac.evlo)
            eqdeps_.append(tr.stats.sac.evdp)
            evids_.append(tr.stats.sac.kevnm)
            if MPI.COMM_WORLD.Get_rank() == 0:
                data_.append(tr.data)
            components_.append(tr.stats.sac.kcmpnm)
        
        theta_phi = [_calthetaphi(stalat, stalon, eqlat, eqlon) 
                     for stalat, stalon, eqlat, eqlon 
                     in zip(lats_, lons_, eqlats_, eqlons_)]
        thetas_ = np.array([x[0] for x in theta_phi], dtype=np.float64)
        phis_ = np.array([x[1] for x in theta_phi], dtype=np.float64)

        dataset_info = pd.DataFrame(dict(
            lats=lats_, lons=lons_, names=names_,
            nets=nets_, thetas=thetas_, phis=phis_,
            eqlats=eqlats_, eqlons=eqlons_,
            eqdeps=eqdeps_, evids=evids_, indices=indices_
        ))
        # drop dupplicate sac files with identical source/receiver
        # values, due to multiple seismic components
        n_before = len(traces)
        dataset_info.drop_duplicates(
            subset=['names', 'nets', 'evids'],
            inplace=True)
        n_after = len(dataset_info)
        if verbose >= 1:
            print('Dropped {} sac files'.format(n_before - n_after))
            
        dataset_info.sort_values(by='evids', inplace=True)
        
        dataset_info.index = list(range(n_after))

        nr = len(dataset_info)
        nrs = dataset_info.groupby('evids').count().lats.values
        dataset_event_info = dataset_info.drop_duplicates(['evids'])
        evids = dataset_event_info.evids.values
        eqlats = dataset_event_info.eqlats.values.astype(np.float64)
        eqlons = dataset_event_info.eqlons.values.astype(np.float64)
        eqdeps = dataset_event_info.eqdeps.values.astype(np.float64)
        r0s = 6371. - eqdeps

        # read event catalog
        # TODO case when event_id is not in the catalog
        cat = read_catalog()
        events_ = cat[np.isin(cat, evids)]
        if len(events_) != len(evids):
            raise RuntimeError('Some events not in the catalog')
        mts = np.array([e.mt for e in events_])
        source_time_functions = np.array(
            [e.source_time_function for e  in events_])
        centroid_times = np.array([e.centroid_time for e in events_])

        events = np.array([
            Event(id, lat, lon, depth, mt, ctime, source_time_function)
            for id, lat, lon, depth, mt, ctime, source_time_function
            in zip(
                evids, eqlats, eqlons, eqdeps,
                mts, centroid_times, source_time_functions)])
        stations = dataset_info.apply(
            lambda x: Station(x.names, x.nets, x.lats, x.lons),
            axis=1).values

        phis = dataset_info.phis.values
        thetas = dataset_info.thetas.values
        lons = dataset_info.lons.values
        lats = dataset_info.lats.values

        if MPI.COMM_WORLD.Get_rank() == 0:
            npts = np.array([len(d) for d in data_], dtype=int).max()
            data_arr = np.zeros((1, 3, nr, npts), dtype=DATA_FLOAT_PREC)
            for ista in range(len(dataset_info.indices.values)):
                component = components_[dataset_info.indices.values[ista]]
                icomp = Component.parse_component(component).value
                try:
                    data_arr[0, icomp, ista] = data_[
                        dataset_info.indices.values[ista]]
                except:
                    n_tmp = len(data_[dataset_info.indices.values[ista]])
                    if n_tmp < npts:
                        tmp_data = np.pad(
                            data_[dataset_info.indices.values[ista]],
                            (0,npts-n_tmp), mode='constant',
                            constant_values=(0,0))
                        data_arr[0, icomp, ista] = tmp_data
                    else:
                        data_arr[0, icomp, ista] = (
                            data_[dataset_info.indices.values[ista]][:npts])

            remaining_traces_indices = (set(indices_)
                - set(dataset_info.indices.values))

            for i in remaining_traces_indices:
                dataset_filt = dataset_info[
                    (dataset_info.evids == evids_[i])
                    & (dataset_info.names == names_[i])
                    & (dataset_info.nets == nets_[i])]
                j = dataset_filt.index.values[0]
                component = components_[i]
                icomp = Component.parse_component(component).value
                try:
                    data_arr[0, icomp, j] = data_[i]
                except:
                    n_tmp = len(data_[i])
                    if n_tmp < npts:
                        tmp_data = np.pad(
                            data_[i],
                            (0,npts-n_tmp), mode='constant',
                            constant_values=(0,0))
                        data_arr[0, icomp, j] = tmp_data
                    else:
                        data_arr[0, icomp, j] = (
                            data_[i][:npts])
        else:
            data_arr = None

        return cls(
            lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events, data_arr, sampling_hz)

    @classmethod
    def dataset_from_sac_process(
            cls, sac_files, windows,
            freq, freq2, filter_type='bandpass',
            verbose=0):
        """Creates a dataset from a list of sac files.
        Data are read from sac files, cut using the time windows,
        and stored in self.data. The sac file data are read and cut
        event by event, which allows to read large dataset.

        This method should be used instead of dataset_from_sac()
        when large amount of data is to be read. It has the same effect
        of using dataset_from_sac() followed by apply_windows(), but
        is much more memory efficient. For instance,10,000 3-components
        records with 20 Hz sampling and 1 hour of
        recording take approx. 138 Gb in memory. The same dataset cut
        in 100 s windows around a single phase (e.g., ScS) takes
        approx 1.9 Gb in memory.

        Args:
            sac_files (list of str): list of paths to sac files.
            windows (list of Window): time windows
            freq (float): minimum filter frequency
            freq2 (float): maximum filter frequency
            filter_type (str): 'bandpass' or 'lowpass'
                (default is 'bandpass')
            verbose (int): 0: quiet, 1: debug.

        Returns:
            Dataset: dataset

        """

        traces = [read(sac_file, headonly=True)[0]
                  for sac_file in sac_files]

        sac_files_by_event = defaultdict(list)
        for sac_file, trace in zip(sac_files, traces):
            sac_files_by_event[trace.stats.sac.kevnm].append(sac_file)

        for i, (event_id, files) in enumerate(sac_files_by_event.items()):
            if i == 0:
                ds = Dataset.dataset_from_sac(
                    files, verbose=verbose, headonly=False)
                ds.filter(freq, freq2, filter_type)
                n_phases = len(set(
                    [(w.phase_name, w.component) for w in windows]))
                npts_max = int(
                        max([w.get_length() for w in windows]) *
                        ds.sampling_hz)
                ds.apply_windows(windows, n_phases, npts_max)
            else:
                ds_tmp = Dataset.dataset_from_sac(
                        files, verbose=verbose, headonly=False)
                ds_tmp.filter(freq, freq2, filter_type)
                n_phases = len(set(
                    [(w.phase_name, w.component) for w in windows]))
                npts_max = int(
                    max([w.get_length() for w in windows]) *
                    ds.sampling_hz)
                ds_tmp.apply_windows(windows, n_phases, npts_max)
                ds.append(ds_tmp)

        return ds

    def append(self, dataset):
        """Append dataset to self.
        """
        assert self.sampling_hz == dataset.sampling_hz
        assert set([event.event_id for event in self.events]).isdisjoint(
            set([event.event_id for event in dataset.events]))

        self.lats = np.hstack([self.lats, dataset.lats])
        self.lons = np.hstack([self.lons, dataset.lons])
        self.phis = np.hstack([self.phis, dataset.phis])
        self.thetas = np.hstack([self.thetas, dataset.thetas])
        self.eqlats = np.hstack([self.eqlats, dataset.eqlats])
        self.eqlons = np.hstack([self.eqlons, dataset.eqlons])
        self.r0s = np.hstack([self.r0s, dataset.r0s])
        self.mts = np.hstack([self.mts, dataset.mts])
        self.nrs = np.hstack([self.nrs, dataset.nrs])
        self.nr = len(self.lats) + len(dataset.lats)
        self.stations = np.hstack([self.stations, dataset.stations])
        self.events = np.hstack([self.events, dataset.events])
        self.data = np.concatenate((self.data, dataset.data), axis=2)
        # self.sampling_hz = sampling_hz
        # self.is_cut = is_cut

    def apply_windows(
            self, windows, n_phase, npts_max, buffer=0.,
            t_before_noise=100., inplace=True):
        """Cut the data using provided windows.

        Args:
            windows (list of Window): time windows.
            n_phase (int): number of distinct seismic phase-component
                pairs: if ScS (SH) and ScS (SV), then n_phase=2.
            npts_max (int): number of time points in the longest window.
            buffer (float): default is 0.
            t_before_noise (float): default is 50.
            inplace (bool): if True, performs the operation in-place
                (i.e., modifies self.data)

        Returns:
            Dataset: if inplace is True, else None.

        """
        if self.data is None:
            return None
        if not inplace:
            ds = self.copy()
        else:
            ds = self
        npts_buffer = int(buffer * ds.sampling_hz)
        data_cut = np.zeros(
            (n_phase, 3, ds.nr, npts_max+2*npts_buffer),
            dtype=DATA_FLOAT_PREC)
        ds.ts_start_end = np.zeros((n_phase, ds.nr, 2), dtype=DATA_FLOAT_PREC)
        ds.noise = np.zeros((n_phase, 3, ds.nr), dtype=DATA_FLOAT_PREC)
        for iev, event in enumerate(ds.events):
            start, end = ds.get_bounds_from_event_index(iev)
            for ista in range(start, end):
                windows_filt = [
                    w for w in windows
                    if w.station == ds.stations[ista]
                    and w.event == event]

                # TODO do it. Current state leads to data dupplication
                # but be careful to add procedure to match the windows

                # Remove windows where only the component differ
                # windows_filt_onecomponent = []
                # for w in windows_filt:
                #     found = False
                #     for w1 in windows_filt_onecomponent:
                #         if (w.event==w1.event
                #             and w.station==w1.station
                #             and w.phase_name == w1.phase_name):
                #             found = True
                #             break
                #     if not found:
                #         windows_filt_onecomponent.append(w)
                # windows_filt = windows_filt_onecomponent

                for iwin, window in enumerate(windows_filt):
                    window_arr = window.to_array()
                    i_start = int((window_arr[0] - buffer) * ds.sampling_hz)
                    i_end = int((window_arr[1] + buffer) * ds.sampling_hz)
                    data_cut[iwin, :, ista, :] = ds.data[0, :,
                        ista, i_start:i_end]
                    ds.ts_start_end[iwin, ista] = window_arr
                    # compute noise
                    i_end_noise = (i_start
                        - int(t_before_noise*ds.sampling_hz))
                    i_start_noise = (i_end_noise - (i_end-i_start))
                    noise_tr = ds.data[0, :, ista, i_start_noise:i_end_noise]
                    noise = np.sqrt(
                        np.sum(noise_tr**2, axis=1) / noise_tr.shape[1])
                    ds.noise[iwin, :, ista] = noise
        ds.data = data_cut
        ds.is_cut = True

        if not inplace:
            return ds
        else:
            return None

    def get_chunks_station(self, n_cores, verbose=0):
        chunk_size = self.nr // n_cores
        dividers = self.nrs / chunk_size
        dividers = Dataset._round_dividers(dividers, n_cores)
        assert np.sum(dividers) == n_cores
        if (dividers == 0).sum() > 0:
            if verbose == 2:
                print('dividers: {}'.format(dividers))
            raise RuntimeError(
                'n_cores ({}) must be >= number of eqs ({})'
                .format(n_cores, len(self.events)))
        counts = self.nrs / dividers
        counts_ = []
        for i in range(len(dividers)):
            counts_.append(Dataset._split(self.nrs[i], counts[i]))
        counts = np.concatenate(counts_)
        assert len(counts.flatten()) == n_cores
        displacements = np.concatenate([[0], counts.cumsum()[:-1]])
        return counts, displacements

    def get_chunks_eq(self, n_cores):
        chunk_size = self.nr // n_cores
        dividers = self.nrs / chunk_size
        dividers = Dataset._round_dividers(dividers, n_cores)
        assert np.sum(dividers) == n_cores
        counts = np.ones(dividers.sum())
        displacements_ = []
        for i, divider in enumerate(dividers):
            disp = np.empty(divider)
            disp.fill(i)
            displacements_.append(disp)
        displacements = np.concatenate(displacements_)
        return counts, displacements

    def get_chunks_mt(self, n_cores):
        counts, displacements = self.get_chunks_eq(n_cores)
        return 9*counts, 9*displacements

    def filter(
            self, freq, freq2=0., type='bandpass', zerophase=False,
            inplace=True):
        """Filter waveforms using obspy.signal.filter.

        Args:
            freq (float): filter frequency.
            freq2 (float): filter maximum frequency.
                For bandpass filters only.
            type (str): type of filter. 'lowpass' or 'bandpass'.
            zerophase (bool): use zero phase filter.
            inplace (bool): if True, performs the operation in-place
                (i.e., modifies self.data).

        Returns:
            Dataset: if inplace is True, else None

        """
        if self.data is None:
            return None
        if type == 'bandpass':
            assert freq2 > freq
        if type not in {'bandpass', 'lowpass'}:
            raise ValueError('Expect "bandpass" or "lowpass" for arg "type"')
        if self.data.shape[3] == 0:
            print('Do nothing')
            return

        if not inplace:
            ds = self.copy()
        else:
            ds = self

        for iwin in range(ds.data.shape[0]):
            for icomp in range(3):
                for ista in range(ds.data.shape[2]):
                    if type == 'lowpass':
                        ds.data[iwin, icomp, ista] = (
                            obspy.signal.filter.lowpass(
                                ds.data[iwin, icomp, ista], freq,
                                df=ds.sampling_hz, zerophase=zerophase))
                    elif type == 'bandpass':
                        ds.data[iwin, icomp, ista] = (
                            obspy.signal.filter.bandpass(
                                ds.data[iwin, icomp, ista], freq, freq2,
                                df=ds.sampling_hz, zerophase=zerophase))
        if not inplace:
            return ds
        else:
            return None

    def get_bounds_from_event_index(self, ievent: int) -> (int, int):
        """Return the start, end indices to slice
        self.stations[start:end].

        Args:
            ievent (int): index of the event as in self.events

        Returns:
            int: index of the first station recording event ievent
            int: index of the last station
        """
        start = self.nrs[:ievent].sum()
        end = start + self.nrs[ievent]
        return start, end

    def set_source_time_functions(self, type, catalog_path=None):
        """Set the catalog for source time functions.
        By default, source time functions specified in the GCMT catalog
        are used.

        Args:
            type (str): 'scardec' or 'user'
            catalog_path: path to a custom catalog.
                Must be specified if type='user'
        """
        if type == 'scardec':
            stf_catalog = STFCatalog.read_scardec()
        elif type == 'user':
            if catalog_name is None:
                raise ValueError(
                    'Expect a path to a catalog when type="user"')
            stf_catalog = STFCatalog.read_from_file(catalog_path)
        else:
            raise ValueError('Expect "scardec" or "user" for arg "type"')
        for i, event in enumerate(self.events):
            if event.event_id in stf_catalog:
                stf = stf_catalog[event.event_id]
                self.events[i].source_time_function = stf

    def plot_event(
            self, ievent, windows=None, align_zero=False,
            component=Component.T, ax=None,
            dist_min=0, dist_max=360, **kwargs):
        """Plot a record section for event ievent.

        Args:
            ievent (int): index of the event as in self.events
            windows (list of Window): time windows used to cut the
                waveforms if specified (default is None)
            align_zero (bool): if True, set the start of time windows
                as t=0 (default is False)
            component (Component): seismic component
                (default is Component.T)
            ax (Axes): matplotlib Axes object
            dist_min (float): minimum epicentral distance (default is 0)
            dist_max (float) maximum epicentral distances (default is 360)
            **kwargs: key-value arguments for the pyplot.plot function

        Returns:
            Figure: matplotlib Figure object
            Axes: matplotlib Axes object

        """
        start, end = self.get_bounds_from_event_index(ievent)
        if ax == None:
            fig, ax = plt.subplots(1)
        else:
            fig = None
        for i in range(start, end):
            distance = self.events[ievent].get_epicentral_distance(
                self.stations[i])
            if distance < dist_min or distance > dist_max:
                continue

            # select corresponding window
            if (windows is not None) and (not self.is_cut):
                windows_tmp = list(filter(
                    lambda w: ((w.station == self.stations[i])
                                and (w.event == self.events[ievent])
                                and (w.component == component)),
                    windows))
                if not windows_tmp:
                    continue
                window = windows_tmp[0].to_array()
                i0 = int(window[0] * self.sampling_hz)
                i1 = int(window[1] * self.sampling_hz)

                data = self.data[0, component.value, i, i0:i1]
                ts = np.linspace(window[0], window[1], len(data))
            else:
                data = self.data[0, component.value, i]
                ts = np.linspace(0, len(data)/self.sampling_hz, len(data))
            if align_zero:
                ts = np.linspace(
                    0, len(data)/self.sampling_hz, len(data))

            norm = np.abs(data).max() * 2.
            norm = norm if norm > 0 else 1.
            ax.plot(ts, data/norm+distance, **kwargs)
            ax.set(xlabel='Time (s)', ylabel='Distance (deg)')
        return fig, ax
    
    @staticmethod
    def _round_dividers(dividers, n_cores):
        dividers_rounded = np.round(dividers).astype(np.int64)
        if np.sum(dividers_rounded) < n_cores:
            dividers_fixed = np.where(
                dividers_rounded == 0, 1, dividers_rounded)
            if np.sum(dividers_fixed) <= n_cores:
                dividers_rounded = dividers_fixed
        dividers_sorted = np.sort(dividers)
        i = -1
        while np.sum(dividers_rounded) != n_cores:
            index_current_max = np.argwhere(dividers == dividers_sorted[i])[0]
            # to avoid problems with identical values
            dividers[index_current_max] = 0.
            dividers_rounded[index_current_max] += 1
            i -= 1
        if (dividers_rounded.sum() == n_cores
            and (dividers_rounded==0).sum() > 0):
            indices_zero = np.argwhere(dividers_rounded == 0)
            for i0 in indices_zero:
                imax = np.argmax(dividers_rounded)
                dividers_rounded[imax] -= 1
                dividers_rounded[i0] += 1
        return dividers_rounded

    @staticmethod
    def _split(size, chunk_size):
        n = size / chunk_size
        n = int(n)
        splits = np.empty(n, dtype=np.int64)
        splits.fill(int(chunk_size))
        splits[-1] = size - (n-1) * int(chunk_size)
        return splits

def filter_sac_files(sac_files, f):
    """Filter sac files using the boolean function f.

    Args:
        sac_files (list of str): paths to sac files
        f (function): (event_id: str, station: Station) -> bool

    Returns:
        list of str: filtered list of paths to sac files

    """
    traces = [read(sac_file, headonly=True)[0]
              for sac_file in sac_files]

    x = lambda tr: (get_event_id(tr), get_station(tr))
    mask = list(map(f, list(map(x, traces))))

    sac_files_filt = [sac_file for i, sac_file in enumerate(sac_files)
                      if mask[i]]
    return sac_files_filt


def filter_abnormal_data(sac_files, f, threshold=5):
    """Filter sac data using the boolean function f.

    Args:
        sac_files (list of str): paths to sac files
        f (function): (event_id: str, station: Station) -> bool
        threshold (float): number of standard deviations of
            the distribution of the log of max of data within which to
            keep the data (default is 5).

    Returns:
        list of str: filtered list of paths to sac files

    """
    traces = [read(sac_file, headonly=False)[0]
              for sac_file in sac_files]
    mask = np.ones(len(traces), dtype='bool')
    maxs = np.zeros(len(traces))

    for i, trace in enumerate(traces):
        if np.isnan(trace.data).any():
            mask[i] = False
        if (trace.data == 0).all():
            mask[i] = False
        maxs[i] = np.max(np.abs(trace.data))

    maxs_filt = maxs[mask]
    log_maxs = np.log(maxs_filt)
    mean_log_max = log_maxs.mean()
    std_log_max = log_maxs.std()

    up = mean_log_max + threshold * std_log_max
    lo = mean_log_max - threshold * std_log_max

    for i, trace in enumerate(traces):
        if mask[i]:
            log_max = np.log(maxs[i])
            if (log_max > up) or (log_max < lo):
                mask[i] = False

    return [file for i, file in enumerate(sac_files) if mask[i]]


def get_event_id(trace):
    """Return event GCMT ID from obspy Trace.

    Args:
        trace (Trace): obspy Trace object

    Returns:
        str: GCMT ID

    """
    evid = trace.stats.sac.kevnm
    return evid

def get_station(trace):
    """Return Station object from obspy Trace.

    Args:
        trace (Trace): obspy Trace object

    Returns:
        Station: station

    """
    sta_nm = trace.stats.sac.kstnm
    sta_net = trace.stats.sac.knetwk
    sta_la = trace.stats.sac.stla
    sta_lo = trace.stats.sac.stlo
    return Station(sta_nm, sta_net, sta_la, sta_lo)
