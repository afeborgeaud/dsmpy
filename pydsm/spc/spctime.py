import numpy as np
import matplotlib.pyplot as plt
import warnings

class SpcTime:
    def __init__(self, tlen, nspc, sampling_hz, omegai,
                 source_time_function=None):
        self.tlen = tlen
        self.nspc = nspc
        self.sampling_hz = sampling_hz
        self.npts = self.find_npts()
        self.lsmooth = self.find_lsmooth()
        self.omegai = omegai
        self.ncomp = 3
        self.source_time_function = source_time_function

    def find_npts(self):
        npts = int(self.tlen * self.sampling_hz)
        pow2 = self.setBitNumber(npts)
        pow2 = pow2 * 2 if pow2 < npts else pow2
        # return pow2
        return npts

    def find_lsmooth(self):
        nspc = self.setBitNumber(self.nspc)
        if nspc < self.nspc:
            nspc *= 2
        lsmooth = self.npts // (nspc * 2)
        i = self.setBitNumber(lsmooth)
        if i < lsmooth:
            i *= 2
        return i

    def setBitNumber(self, n):
        n |= n >> 1
        n |= n >> 2
        n |= n >> 4
        n |= n >> 8
        n |= n >> 16
        n = n + 1
        return (n >> 1)

    def to_time_domain(self, spc):
        assert self.npts % 2 == 0
        nnp = self.npts // 2
        uspc = np.pad(spc, pad_width=(0, nnp - len(spc)), mode='constant',
                      constant_values=0)
        uspc_conj = np.pad(np.flip(spc[1:]).conjugate(),
                           pad_width=(nnp - len(spc) + 1, 0), mode='constant',
                           constant_values=0)
        uspc = np.concatenate((uspc, uspc_conj))

        ureal = np.real(np.fft.ifft(uspc))
        return ureal.astype(np.float64)

    def apply_growing_exponential(self, u):
        c = self.omegai * self.tlen
        x = np.linspace(0, c, self.npts, endpoint=False)
        a = np.exp(x)
        np.multiply(u, a, out=u)

    def apply_amplitude_correction(self, u):
        c = self.npts * 1e3 / self.tlen
        u *= c

    def convolve(self, spc):
        if self.source_time_function is not None:
            stf_freq = (self.source_time_function.
                get_source_time_function_frequency_domain(
                    self.tlen, self.nspc))
            return spc * stf_freq
        else:
            return spc

    def spctime(self, spcs):
        ''' spcs is the output of pyDSM.
            spcs.shape = (3, nr, imax)
            return:
            u: time-domain version of spcs
        '''
        nr = spcs.shape[1]
        u = np.zeros((3, nr, self.npts))
        for icomp in range(3):
            for ir in range(nr):
                spcs_conv = self.convolve(spcs[icomp, ir])
                u[icomp, ir, :] = self.to_time_domain(spcs_conv)
                self.apply_growing_exponential(u[icomp, ir])
                self.apply_amplitude_correction(u[icomp, ir])
        return u
