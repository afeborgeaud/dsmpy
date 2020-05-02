import numpy as np
import matplotlib.pyplot as plt

class SpcTime:
    
    def __init__(self, tlen, nspc, sampling_hz, omegai):
        self.tlen = tlen
        self.nspc = nspc
        self.sampling_hz = sampling_hz
        self.npts = self.find_npts()
        self.lsmooth = self.find_lsmooth()
        self.omegai = omegai

    def find_npts(self):
        npts = int(self.tlen * self.sampling_hz)
        pow2 = self.setBitNumber(npts)
        pow2 = pow2 * 2 if pow2 < npts else pow2
        return pow2

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
        n |= n>>1
        n |= n>>2
        n |= n>>4  
        n |= n>>8
        n |= n>>16
        n = n + 1
        return (n >> 1) 
    
    def to_time_domain(self, spc):
        nnp = self.npts // 2
        uspc = np.pad(spc, pad_width=(0, nnp-len(spc)), mode='constant',
            constant_values=0)
        uspc_conj = np.pad(np.conjugate(np.flip(spc)), pad_width=(nnp-len(spc), 0), mode='constant',
            constant_values=0)
        uspc = np.concatenate((uspc, uspc_conj))

        ureal = np.real(np.fft.ifft(uspc))[:]

        return ureal.astype(np.float64)

    def apply_growing_exponential(self, u):
        c = self.omegai * self.tlen / self.npts
        x = np.linspace(0, c, self.npts)
        a = np.exp(x)
        return u * a

    def apply_amplitude_correction(self, u):
        c = self.npts * 1e3 / self.tlen
        return u * c

    def spctime(self, spc):
        u = self.to_time_domain(spc)
        u = self.apply_growing_exponential(u)
        u = self.apply_amplitude_correction(u)
        return u