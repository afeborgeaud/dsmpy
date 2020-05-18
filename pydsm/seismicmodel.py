import numpy as np
from pydsm._tish import parameters

class SeismicModel:
    """Represent a seismic Earth model for computation using DSM.

    Args:
        vrmin (ndarray): 1D array containing
            lower bounds of layered structure [km]
        vrmax (ndarray): 1D array containing
            upper bounds of layered structure [km]
        rho (ndarray): 2D array specifying density using
            3-degree polynomials for each layer [g/cm^3]
        vpv (ndarray): 2D array specifying V_PV using
            3-degree polynomials for each layer [km/s]
        vph (ndarray): V_PH [km/s]
        vsh (ndarray): V_SH [km/s]
        eta (ndarray): radial anisotropy []
        qmu (ndarray): shear anelasic factor []
        qkappa (ndarray): bulk anelasic factor []
    
    Attributes:
    """

    def __init__(
            self, vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa):
        self._vrmin = vrmin
        self._vrmax = vrmax
        self._rho = rho
        self._vpv = vpv
        self._vph = vph
        self._vsv = vsv
        self._vsh = vsh
        self._eta = eta
        self._qmu = qmu
        self._qkappa = qkappa
        self._nzone = rho.shape[1]

    def get_rho(self):
        return np.pad(
            self._rho, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_vpv(self):
        return np.pad(
            self._vpv, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_vph(self):
        return np.pad(
            self._vph, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_vsv(self):
        return np.pad(
            self._vsv, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_vsh(self):
        return np.pad(
            self._vsh, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_eta(self):
        return np.pad(
            self._eta, ((0, 0), (0, parameters['maxnzone']-self._nzone)),
            mode='constant', constant_values=0)
    def get_vrmin(self):
        return np.pad(
            self._vrmin, (0, parameters['maxnzone']-self._nzone),
            mode='constant', constant_values=0)
    def get_vrmax(self):
        return np.pad(
            self._vrmax, (0, parameters['maxnzone']-self._nzone),
            mode='constant', constant_values=0)
    def get_qmu(self):
        return np.pad(
            self._qmu, (0, parameters['maxnzone']-self._nzone),
            mode='constant', constant_values=0)
    def get_qkappa(self):
        return np.pad(
            self._qkappa, (0, parameters['maxnzone']-self._nzone),
            mode='constant', constant_values=0)

    @classmethod
    def ak135(cls):
        """Return model AK135.

        References:
            Kennett et al. (1995)
        """
        vrmin = np.array([
            0, 1217.5, 3479.5, 3631, 5611, 5711,
            5961, 6161, 6251, 6336.6, 6351
        ], dtype=np.float64)
        vrmax = np.array([
            1217.5, 3479.5, 3631, 5611, 5711,
            5961, 6161, 6251, 6336.6, 6351, 6371
        ], dtype=np.float64)
        rho = np.array([
            [13.0885, 0, -8.8381, 0],
            [12.5815, -1.2638, -3.6426, -5.5281],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [5.3197, -1.4836, 0, 0],
            [11.2494, -8.0298, 0, 0],
            [7.1089, -3.8045, 0, 0],
            [2.691, 0.6924, 0, 0],
            [2.691, 0.6924, 0, 0],
            [2.9, 0, 0, 0],
        ], dtype=np.float64).T
        vpv = np.array([
            [11.261692, 0.028794, -6.627846, 0],
            [10.118851, 3.457774, -13.434875, 0],
            [13.908244, -0.45417, 0, 0],
            [24.138794, -37.097655, 46.631994, -24.272115],
            [25.969838, -16.934118, 0, 0],
            [29.38896, -21.40656, 0, 0],
            [30.78765, -23.25415, 0, 0],
            [25.413889, -17.697222, 0, 0],
            [8.785412, -0.749529, 0, 0,],
            [6.5, 0.0, 0.0, 0.0],
            [5.8, 0.0, 0.0, 0.0]
        ], dtype=np.float64).T
        vph = np.array([
            [11.261692, 0.028794, -6.627846, 0],
            [10.118851, 3.457774, -13.434875, 0],
            [13.908244, -0.45417, 0, 0],
            [24.138794, -37.097655, 46.631994, -24.272115],
            [25.969838, -16.934118, 0, 0],
            [29.38896, -21.40656, 0, 0],
            [30.78765, -23.25415, 0, 0],
            [25.413889, -17.697222, 0, 0],
            [8.785412, -0.749529, 0, 0,],
            [6.5, 0, 0, 0],
            [5.8, 0, 0, 0]
        ], dtype=np.float64).T
        vsv = np.array([
            [3.667865, -0.001345, -4.440915, 0],
            [0, 0, 0, 0],
            [8.018341, -1.349895, 0, 0],
            [12.213901, -18.573085, 24.557329, -12.728015],
            [20.208945, -15.895645, 0, 0],
            [17.71732, -13.50652, 0, 0],
            [15.212335, -11.053685, 0, 0],
            [5.7502, -1.2742, 0, 0],
            [5.970824, -1.499059, 0, 0],
            [3.85, 0, 0, 0],
            [3.46, 0, 0, 0]
        ], dtype=np.float64).T
        vsh = np.array([
            [3.667865, -0.001345, -4.440915, 0],
            [0, 0, 0, 0],
            [8.018341, -1.349895, 0, 0],
            [12.213901, -18.573085, 24.557329, -12.728015],
            [20.208945, -15.895645, 0, 0],
            [17.71732, -13.50652, 0, 0],
            [15.212335, -11.053685, 0, 0],
            [5.7502, -1.2742, 0, 0],
            [5.970824, -1.499059, 0, 0],
            [3.85, 0, 0, 0],
            [3.46, 0, 0, 0]
        ], dtype=np.float64).T
        eta = np.array([
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0]
        ], dtype=np.float64).T
        qmu = np.array([
            84.6, -1, 312, 312, 312, 143,
            143, 80, 600, 600, 600
        ], dtype=np.float64)
        qkappa = np.array([
            1327.7, 57823, 57823, 57823, 57823, 57823,
            57823, 57823, 57823, 57823, 57823
        ], dtype=np.float64)
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)

    @classmethod
    def prem(cls):
        """Return the Preliminary Reference Earth Model (PREM).

        References:
            Dziewonski and Anderson (1981)
        """
        vrmin = np.array([
            0, 1221.5, 3480, 3630, 5600, 5701, 5771,
            5971, 6151, 6291, 6346.6, 6356
            ], dtype=np.float64)
        vrmax = np.array([
            1221.5, 3480, 3630, 5600, 5701, 5771,
            5971, 6151, 6291, 6346.6, 6356, 6371
            ], dtype=np.float64)
        rho = np.array([
            [13.0885, 0, -8.8381, 0],
            [12.5815, -1.2638, -3.6426, -5.5281],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [7.9565, -6.4761, 5.5283, -3.0807],
            [5.3197, -1.4836, 0, 0],
            [11.2494, -8.0298, 0, 0],
            [7.1089, -3.8045, 0, 0],
            [2.691, 0.6924, 0, 0],
            [2.691, 0.6924, 0, 0],
            [2.9, 0, 0, 0],
            [2.6, 0, 0, 0]
        ], dtype=np.float64).T
        vpv = np.array([
            [11.2622, 0, -6.364, 0],
            [11.0487, -4.0362, 4.8023, -13.5732],
            [15.3891, -5.3181, 5.5242, -2.5514],
            [24.952, -40.4673, 51.4832, -26.6419],
            [29.2766, -23.6027, 5.5242, -2.5514],
            [19.0957, -9.8672, 0, 0],
            [39.7027, -32.6166, 0, 0],
            [20.3926, -12.2569, 0, 0],
            [0.8317, 7.218, 0, 0],
            [0.8317, 7.218, 0, 0],
            [6.8, 0, 0, 0],
            [5.8, 0, 0, 0]
        ], dtype=np.float64).T
        vph = np.array([
            [11.2622, 0, -6.364, 0],
             [11.0487, -4.0362, 4.8023, -13.5732],
            [15.3891, -5.3181, 5.5242, -2.5514],
             [24.952, -40.4673, 51.4832, -26.6419],
            [29.2766, -23.6027, 5.5242, -2.5514],
             [19.0957, -9.8672, 0, 0],
              [39.7027, -32.6166, 0, 0],
            [20.3926, -12.2569, 0, 0],
             [3.5908, 4.6172, 0, 0],
              [3.5908, 4.6172, 0, 0],
               [6.8, 0, 0, 0],
            [5.8, 0, 0, 0]
        ], dtype=np.float64).T
        vsv = np.array([
            [3.6678, 0, -4.4475, 0],
            [0, 0, 0, 0],
            [6.9254, 1.4672, -2.0834, 0.9783],
            [11.1671, -13.7818, 17.4575, -9.2777],
            [22.3459, -17.2473, -2.0834, 0.9783],
            [9.9839, -4.9324, 0, 0],
            [22.3512, -18.5856, 0, 0],
            [8.9496, -4.4597, 0, 0],
            [5.8582, -1.4678, 0, 0], [5.8582, -1.4678, 0, 0], [3.9, 0, 0, 0],
            [3.2, 0, 0, 0]
        ], dtype=np.float64).T
        vsh = np.array([
            [3.6678, 0, -4.4475, 0],
            [0, 0, 0, 0],
            [6.9254, 1.4672, -2.0834, 0.9783],
            [11.1671, -13.7818, 17.4575, -9.2777],
            [22.3459, -17.2473, -2.0834, 0.9783],
            [9.9839, -4.9324, 0, 0],
            [22.3512, -18.5856, 0, 0],
            [8.9496, -4.4597, 0, 0],
            [-1.0839, 5.7176, 0, 0],
            [-1.0839, 5.7176, 0, 0],
            [3.9, 0, 0, 0],
            [3.2, 0, 0, 0]
        ], dtype=np.float64).T
        eta = np.array([
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [3.3687, -2.4778, 0, 0],
            [3.3687, -2.4778, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0]
        ], dtype=np.float64).T
        qmu = np.array([
            84.6, 1e12, 312, 312, 312, 143,
            143, 143, 80, 600, 600, 600
        ], dtype=np.float64)
        qkappa = np.array([
            1327.7, 57823, 57823, 57823, 57823,
            57823, 57823, 57823, 57823, 57823, 57823, 57823
        ], dtype=np.float64)
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)

    @classmethod
    def iasp91(cls): # TODO
        """Return model IAS91.

        References:
        """
        vrmin = None
        vrmax = None
        rho = None
        vpv = None
        vph = None
        vsv = None
        vsh = None
        eta = None
        qmu = None
        qkappa = None
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)
