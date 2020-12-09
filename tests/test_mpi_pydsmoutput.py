from mpi4py import MPI
from dsmpy.utils.cmtcatalog import read_catalog
from dsmpy.dsm import PyDSMInput, compute
from dsmpy.seismicmodel import SeismicModel
from dsmpy.station import Station
from dsmpy.event import Event

def compute_output(tlen=1638.4, nspc=64, sampling_hz=20, mode=2):
    catalog = read_catalog()
    event = Event.event_from_catalog(
        catalog, '200707211534A')
    stations = [
        Station(
            '{:03d}'.format(i), 'DSM', event.latitude, event.longitude+i)
        for i in range(12,36)]
    
    model = SeismicModel.ak135()
    pydsm_input = PyDSMInput.input_from_arrays(
        event, stations, model, tlen, nspc, sampling_hz)
    pydsm_output = compute(pydsm_input, mode=mode)
    
    return pydsm_output

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        output = compute_output()
        outputs = [[output for i in range(2)] for i in range(2)]
    else:
        outputs = None

    outputs = comm.scatter(outputs, root=0)
    print('{} {}'.format(rank, len(outputs)))
    print('{} {}'.format(rank, outputs[0].event))

    # print('{} {}'.format(rank, output.event))