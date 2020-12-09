import warnings
import numpy as np

class SourceTimeFunction:
    """Represent an earthquake source time function.
    
    Args:
        type (str): 'triangle' or 'box car'
        half_duration (float): half duration of the source time function
    """
    _types = {'triangle', 'box car'}

    def __init__(
            self, type: str, half_duration: float,
            amp_corr=1.):
        self.type = type
        self.half_duration = half_duration
        self.amp_corr = amp_corr

        if self.type not in self._types:
            raise RuntimeError('{} not implemented yet'.format(self.type))

    def get_source_time_function_frequency_domain(
            self, tlen, nspc, amp_corr=1.):
        if self.type == 'triangle':
            return SourceTimeFunction.triangle(
                self.half_duration, tlen, nspc, amp_corr)
        elif self.type == 'box car':
            return SourceTimeFunction.boxcar(
                self.half_duration, tlen, nspc, amp_corr)
        else:
            warnings.warn('{} not implemented yet'.format(self.type))
            return None

    @staticmethod
    def triangle(half_duration, tlen, nspc, amp_corr=1.):
        deltaF = 1 / tlen
        constant = 2 * np.pi * deltaF * half_duration
        stf = np.zeros(nspc + 1, dtype=np.complex128)
        for i in range(nspc):
            omega_tau = (i + 1) * constant
            stf[i] = complex((2 - 2 * np.cos(omega_tau))
                             / (omega_tau * omega_tau))
        stf *= amp_corr
        return stf

    @staticmethod
    def boxcar(half_duration, tlen, nspc, amp_corr=1.):
        deltaF = 1 / tlen
        constant = 2 * np.pi * deltaF * half_duration
        stf = np.zeros(nspc + 1, dtype=np.complex128)
        for i in range(nspc):
            omega_tau = (i + 1) * constant
            stf[i] = complex(np.sin(omega_tau) / omega_tau)
        stf *= amp_corr
        return stf