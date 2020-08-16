from pydsm.dsm import PyDSMInput
import numpy as np
from obspy import read
import obspy.signal.filter
import pandas as pd
from pydsm.event import Event, MomentTensor
from pydsm.station import Station
from pydsm._tish import _calthetaphi
from pydsm import root_resources
from pydsm.utils.cmtcatalog import read_catalog

class Dataset:
    """Represent a dataset of events and stations.
    """
    def __init__(
            self, lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events, data=[], sampling=20):
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
        self.data = data
        self.sampling = sampling

    @classmethod
    def dataset_from_files(cls, parameter_files, mode=1):
        pydsm_inputs = [PyDSMInput.input_from_file(file, mode=mode)
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
        mts = np.concatenate([input.mt for input in pydsm_inputs])
        nrs = np.array([input.nr for input in pydsm_inputs])

        stations = np.concatenate([input.stations
                                        for input in pydsm_inputs])
        events = np.array([input.event
                                for input in pydsm_inputs])

        return cls(
            lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events)
    
    @classmethod
    def dataset_from_arrays(cls, events, stations, sampling_hz):
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
            r0s, mts, nrs, stations_, events_, data=[], sampling=sampling_hz)

    @classmethod
    def dataset_from_sac(cls, sac_files, verbose=0, headonly=True):
        headers = [read(sac_file, headonly=headonly)[0]
                   for sac_file in sac_files]

        sampling = headers[0].stats.sampling_rate

        lats = []
        lons = []
        names = []
        nets = []
        eqlats = []
        eqlons = []
        eqdeps = []
        evids = []
        data = []
        indices = list(range(len(headers)))

        for h in headers:
            lats.append(h.stats.sac.stla)
            lons.append(h.stats.sac.stlo)
            names.append(h.stats.sac.kstnm)
            nets.append(h.stats.sac.knetwk)
            eqlats.append(h.stats.sac.evla)
            eqlons.append(h.stats.sac.evlo)
            eqdeps.append(h.stats.sac.evdp)
            evids.append(h.stats.sac.kevnm)
            tr_filt = h.filter('lowpass', freq=1., zerophase=True)
            data.append(tr_filt.data)
        
        dataset_info = pd.DataFrame(dict(
            lats=lats, lons=lons, names=names,
            nets=nets, eqlats=eqlats, eqlons=eqlons,
            eqdeps=eqdeps, evids=evids, indices=indices
        ))
        # drop dupplicate sac files with identical source/receiver
        # values, due to multiple seismic components
        n_before = len(headers)
        dataset_info.drop_duplicates(inplace=True)
        n_after = len(dataset_info)
        if verbose >= 1:
            print('Dropped {} sac files'.format(n_before - n_after))
            
        dataset_info.sort_values(by='evids', inplace=True)

        theta_phi = [_calthetaphi(stalat, stalon, eqlat, eqlon) 
                     for stalat, stalon, eqlat, eqlon 
                     in zip(lats, lons, eqlats, eqlons)]
        thetas = np.array([x[0] for x in theta_phi], dtype=np.float64)
        phis = np.array([x[1] for x in theta_phi], dtype=np.float64)

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

        events = [
            Event(id, lat, lon, depth, mt, ctime, source_time_function)
            for id, lat, lon, depth, mt, ctime, source_time_function
            in zip(
                evids, eqlats, eqlons, eqdeps,
                mts, centroid_times, source_time_functions)]
        stations = dataset_info.apply(
            lambda x: Station(x.names, x.nets, x.lats, x.lons),
            axis=1).values

        lons = np.array(lons, dtype=np.float64)
        lats = np.array(lats, dtype=np.float64)

        data_ = np.array(
            [data[i] for i in dataset_info.indices.values],
            dtype=np.float64)

        return cls(
            lats, lons, phis, thetas, eqlats, eqlons,
            r0s, mts, nrs, stations, events, data_, sampling)

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
        return 9*counts, displacements

    def filter(self, freq, freq1=0., mode='lowpass', zerophase=False):
        for i in range(len(self.data)):
            if mode == 'lowpass':
                self.data[i] = obspy.signal.filter.lowpass(
                    self.data[i], freq, df=self.sampling, zerophase=zerophase)
    
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
            index_current_max = np.argwhere(dividers == dividers_sorted[i])
            dividers_rounded[index_current_max] += 1
            i -= 1
        return dividers_rounded

    @staticmethod
    def _split(size, chunk_size):
        n = size / chunk_size
        n = int(n)
        splits = np.empty(n, dtype=np.int64)
        splits.fill(int(chunk_size))
        splits[-1] = size - (n-1) * int(chunk_size)
        return splits

    def get_bounds_from_event_index(self, ievent):
        start = self.nrs[:ievent].sum()
        end = start + self.nrs[ievent]
        return start, end

