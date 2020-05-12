

class Dataset:
    """Represent a dataset of events and stations.
    """

    def __init__(self, event_array, station_array):
        self.event_array = event_array
        self.station_array = station_array