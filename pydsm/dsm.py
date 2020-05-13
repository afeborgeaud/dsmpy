import sys
from pydsm._tish import _tish
from pydsm._tipsv import _pinput, _tipsv
from pydsm.spc import spctime
import numpy as np
from mpi4py import MPI
import time

class PyDSMOutput:
    def __init__(self, spcs, dsm_input):
        self.spcs = spcs
        self.input = dsm_input
        self.stations = dsm_input.stations
        self.event = dsm_input.event
        self.components = ('Z', 'R', 'T')
        self.dt = 1 / self.input.sampling_hz

    def to_time_domain(self):
        spct = spctime.SpcTime(self.input.tlen, self.input.nspc,
                               self.input.sampling_hz, self.input.omegai,
                               self.input.source_time_function)
        us = spct.spctime(self.spcs)
        self.us = us
        self.ts = np.linspace(0, self.input.tlen,
                              spct.npts, endpoint=False)

    def get_nr(self):
        return self.input.nr


class DSMInput:
    def __init__(
            self, re, ratc, ratl, tlen, nspc, omegai, imin, imax, nzone,
            vrmin, vrmax, rho, vpv, vph, vsv, vsh, eta, qmu, qkappa,
            r0, eqlat, eqlon, mt, nr, theta, phi, lat, lon, output):
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

    def _get_scalar_dict(self):
        return dict(
            re=self.re, ratc=self.ratc, ratl=self.ratl, tlen=self.tlen,
            nspc=self.nspc, omegai=self.omegai, imin=self.imin,
            imax=self.imax, nzone=self.nzone, r0=self.r0, 
            eqlat=self.eqlat, eqlon=self.eqlon, nr=self.nr,
            max_nzone=len(self.vrmin), max_nr=len(self.lat))

    @classmethod
    def input_from_file(cls, parameter_file):
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

        return cls(
            re, ratc, ratl, tlen, nspc,
            omegai, imin, imax, nzone, vrmin, vrmax,
            rho, vpv, vph, vsv, vsh, eta, qmu, qkappa, r0, eqlat, eqlon,
            mt, nr, theta, phi, lat, lon, output)

    @classmethod
    def input_from_arrays(cls, event, stations,
                          seismicmodel, tlen, nspc):
        pass  # TODO

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
        nzone = self.nzone - 2
        vrmin = np.pad(self.vrmin[2:], (0,2), constant_values=0)
        vrmax = np.pad(self.vrmax[2:], (0,2), constant_values=0)
        qmu = np.pad(self.qmu[2:], (0,2), constant_values=0)
        npad = ((0, 0), (0, 2))
        rho = np.pad(self.rho[:,2:], npad, constant_values=0)
        vsv = np.pad(self.vsv[:,2:], npad, constant_values=0)
        vsh = np.pad(self.vsh[:,2:], npad, constant_values=0)

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
    def __init__(
            self, dsm_input, sampling_hz=None,
            source_time_function=None):
        super().__init__(
            *dsm_input.get_inputs_for_tipsv())
        self.source_time_function = source_time_function
        self.sampling_hz = self.find_optimal_sampling_hz(sampling_hz)
        self.stations = self._parse_stations()
        self.event = self._parse_event()

    @classmethod
    def input_from_file(cls, parameter_file,
                        sampling_hz=None, source_time_function=None):
        dsm_input = DSMInput.input_from_file(parameter_file)
        pydsm_input = cls(dsm_input, sampling_hz,
                                 source_time_function)
        return pydsm_input

    @classmethod
    def input_from_arrays(
            cls, event, stations,
            seismicmodel, tlen, nspc, source_time_function,
            sampling_hz):
        dsm_input = super().input_from_arrays(event, stations,
                                              seismicmodel, tlen, nspc)
        pydsm_input = cls(dsm_input, sampling_hz,
                                 source_time_function)
        return pydsm_input

    def find_optimal_sampling_hz(self, sampling_hz):
        if sampling_hz is not None:
            return sampling_hz
        else:
            Tmin = (self.tlen
                    / self.imax)
            optimal_sampling_hz = 40 / Tmin  # TODO: check this scaling
            return optimal_sampling_hz

    def set_sourcetimefunction(self, sourcetimefunction):
        self.sourcetimefunction = sourcetimefunction

    def _parse_stations(self):
        stations = []
        for i in range(self.nr):
            name, net = self.output[i].tostring(). \
                decode('utf-8').split('/')[-1].split('.')[0].split('_')
            station = Station(name, net, self.lat[i], self.lon[i])
            stations.append(station)
        return tuple(stations)

    def _parse_event(self):
        eventID = self.output[0].tostring(). \
            decode('utf-8').split('/')[-1].split('.')[1]
        if eventID[-2:] == 'SH':
            eventID = eventID[:-2]
        elif eventID[-3:] == 'PSV':
            eventID = eventID[:-3]
        else:
            raise RuntimeError('{}'.format(eventID))
        event = Event(eventID, self.eqlat, self.eqlon,
                      6371. - self.r0, self.mt)
        return event


class Station:
    def __init__(self, name: str, network: str,
            latitude: float, longitude: float):
        self.name = name
        self.network = network
        self.latitude = latitude
        self.longitude = longitude

    def __repr__(self):
        return self.name + '_' + self.network


class Event:
    """Represent an earthquake point-source.

    Args:
        eventID (str): GCMT event name
        latitude (float) centroid geographic latitude 
            [-90, 90] in degree
        longitude (float) centroid longitude 
            [-180, 180] in degree
        depth (float) centroid depth in km
        mt (ndarray(3, 3)) moment tensor
    """

    def __init__(self, eventID, latitude, longitude, depth, mt):
        self.eventID = eventID
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.mt = mt

def compute(pydsm_input, mode=0, write_to_file=False):
    """Compute spectra using DSM.

    Args:
        dsm_input (PyDSMInput): inputs for DSM
        mode (int): computation mode. 0: both, 1: P-SV, 2: SH
        write_to_file (bool): write spetrum files to disk
            as specified in dsm_input.output

    Returns:
        dsm_output (DSMOutput): object containing spectra and
            statations/source information
    
    Note:
        SH and P-SV spectra are summed by default. Using only P-SV
        or SH results in non-physical waves and should be avoided.
        See Kawai et al. (2006) for details.
    """
    print('compute SH')
    sh_spcs = _tish(*pydsm_input.get_inputs_for_tish(),
        write_to_file)
    # FIXME memory error in tipsv.tipsv
    #print('compute PSV')
    #psv_spcs = _tipsv(*pydsm_input.get_inputs_for_tipsv(),
    #    write_to_file)
    #spcs = sh_spcs + psv_spcs
    spcs = sh_spcs
    dsm_output = PyDSMOutput(spcs, pydsm_input)
    return dsm_output

def compute_parallel(
        pydsm_input, mode=0, write_to_file=False):
    """Compute spectra using DSM with data parallelization.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_cores = comm.Get_size()
    
    if rank == 0:
        scalar_dict = pydsm_input._get_scalar_dict()
    else:
        scalar_dict = None
    scalar_dict = comm.bcast(scalar_dict, root=0)

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
        comm.send(chunk_sizes[-1], dest=n_cores-1, tag=11)
    if rank == n_cores-1:
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

    lat_local = np.pad(lat_local, (0, scalar_dict['max_nr']-chunk_size),
                       mode='constant', constant_values=0)
    lon_local = np.pad(lon_local, (0, scalar_dict['max_nr']-chunk_size),
                       mode='constant', constant_values=0)
    phi_local = np.pad(phi_local, (0, scalar_dict['max_nr']-chunk_size),
                       mode='constant', constant_values=0)
    theta_local = np.pad(theta_local, (0, scalar_dict['max_nr']-chunk_size),
                       mode='constant', constant_values=0)

    scalar_dict['nr'] = chunk_size
    input_local = DSMInput.input_from_dict_and_arrays(
        scalar_dict, vrmin, vrmax, rho, vpv, vph,
        vsv, vsh, eta, qmu, qkappa, mt, lat_local,
        lon_local, phi_local, theta_local)

    start_time = time.time()
    spcs_local = _tish(*input_local.get_inputs_for_tish(),
                         write_to_file=False)
    end_time = time.time()
    print('{} paths: processor {} in {} s'.format(input_local.nr, rank, end_time - start_time))
    
    # TODO change the order of outputu in DSM 
    # to have nr as the last dimension
    spcs_local = np.array(spcs_local.transpose(0, 2, 1), order='F')

    if rank == 0:
        spcs_gathered = np.empty((3, (pydsm_input.nspc+1), pydsm_input.nr),
                                dtype=np.complex128, order='F')
    else:
        spcs_gathered = None

    if rank == 0:
        spcs_chunk_sizes = tuple([3 * size * (scalar_dict['nspc']+1)
                              for size in chunk_sizes])
        spcs_start_indices = tuple([3 * i * (scalar_dict['nspc']+1)
                              for i in start_indices])
    else:
        spcs_chunk_sizes = None
        spcs_start_indices = None

    comm.Barrier()
    comm.Gatherv(spcs_local, 
                [spcs_gathered, spcs_chunk_sizes, spcs_start_indices,
                 MPI.DOUBLE_COMPLEX],
                root=0)

    return spcs_gathered, comm

def _get_chunk_start_indices(nr, n_cores):
    chunk_size = nr // n_cores
    start_indices = [i * chunk_size for i in range(n_cores)]
    chunk_sizes = [chunk_size for i in range(n_cores-1)]
    chunk_sizes += [nr - start_indices[-1]]
    return tuple(start_indices), tuple(chunk_sizes)