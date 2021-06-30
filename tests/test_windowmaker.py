from dsmpy.windowmaker import WindowMaker
from dsmpy.window import Window
from dsmpy.event import Event
from dsmpy.station import Station
from dsmpy.component import Component
import matplotlib.pyplot as plt

def test_trim_windows():
    event = Event('event1', 0., 0., 100., None, None, None)
    stations = [
        Station(f'{i:03d}', 'SYN', 0., i)
        for i in range(70, 84)
    ]
    windows = WindowMaker.compute(
        event, stations, 'prem', ['S', 'ScS', 'sS'], [Component.T],
        t_before=20, t_after=40)
    windows_ScS = [w for w in windows if w.phase_name == 'ScS']
    windows_S_sS = [w for w in windows
                      if w.phase_name in {'S', 'sS'}]
    windows_S_sS_trim = WindowMaker.set_limit(
        windows_S_sS, 5., 15., inplace=False)
    trimmed_windows_ScS = WindowMaker.trim_windows(
        windows_ScS, windows_S_sS_trim)

    plt.figure(figsize=(5, 7))
    for i, w in enumerate(windows_ScS):
        d = w.get_epicentral_distance()
        label = 'ScS' if i == 0 else None
        plt.plot(w.to_array(), [d, d], color='black', label=label)
    for i, w in enumerate(trimmed_windows_ScS):
        d = w.get_epicentral_distance() + 0.1
        label = 'ScS trimmed' if i == 0 else None
        plt.plot(w.to_array(), [d, d], color='blue', label=label)
    for i, w in enumerate(windows_S_sS_trim):
        d = w.get_epicentral_distance() - 0.1
        label = 'S or sS' if i == 0 else None
        plt.plot(w.to_array(), [d, d], color='red', label=label)
    plt.xlabel('Time (s)')
    plt.ylabel('Distance (deg)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    test_trim_windows()
