from dsmpy.event import Event
from dsmpy.utils.cmtcatalog import read_catalog


def test_event_from_catalog():
    cat = read_catalog()
    event = Event.event_from_catalog(
        cat, '200707211534A')
    assert event.event_id == '200707211534A'


if __name__ == '__main__':
    test_event_from_catalog()