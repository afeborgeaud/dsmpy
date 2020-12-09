class Station:
    """Represent a seismic station.

    Args:
        name (str): station name
        network (str): network code
        latitude (float): geographic latitude
        longitude (float): geographic longitude

    Attributes:
        name (str): station name
        network (str): network code
        latitude (float): geographic latitude
        longitude (float): geographic longitude
    """

    def __init__(self, name: str, network: str,
                 latitude: float, longitude: float):
        self.name = name
        self.network = network
        self.latitude = latitude
        self.longitude = longitude

    def __repr__(self):
        return self.name + '_' + self.network

    def __eq__(self, other):
        if self.__repr__() == other:
            return True
        else:
            return False