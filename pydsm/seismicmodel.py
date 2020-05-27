import numpy as np
from pydsm._tish import parameters
import bisect
from pydsm.modelparameters import ModelParameters, ParameterType
import matplotlib.pyplot as plt
import pickle

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
    def model_from_name(cls, model_name):
        """Return SeismicModel from its identifier.
        Supported models are:
        - ak135
        - prem
        """
        if model_name == 'ak135':
            return cls.ak135()
        elif model_name == 'prem':
            return cls.prem()
        else:
            raise KeyError('{} not implemented yet'.format(model_name))

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
        raise NotImplementedError("IASP91 not yet implemented")
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

    def boxcar_mesh(self, model_parameters):
        """
        Args:
            nodes (ndarray): nodes of the boxcar mesh
        Returns:
            model (SeismicModel): copy of self with added nodes
            mesh (SeismicModel): mesh with boxcar polynomials
        """
        model = self.__copy__()
        for node in model_parameters.get_nodes():
            model = model._add_boundary(node)
        mesh = model.__copy__()
        for i in range(mesh._nzone):
            mesh._set_all_layers(i, np.zeros(4, dtype=np.float64), 0.)
        nodes = model_parameters.get_nodes()
        for i in range(model_parameters._n_nodes-1):
            indexes = (set(np.where(mesh._vrmin < nodes[i+1])[0])
                & set(np.where(mesh._vrmin >= nodes[i])[0]))
            for index in indexes:
                mesh._set_all_layers(
                    index, np.array([1, 0, 0, 0], dtype=np.float64), 1.)
        model._nzone = model._rho.shape[1]
        mesh._nzone = mesh._rho.shape[1]
        return model, mesh

    def multiply(self, nodes, values):
        assert values.shape == (len(nodes)-1, 8)
        assert values.dtype == np.float64
        mesh = self.__copy__()
        for i in range(len(nodes)-1):
            indexes = (set(np.where(self._vrmin < nodes[i+1])[0])
                & set(np.where(self._vrmin >= nodes[i])[0]))
            for index in indexes:
                mesh._rho[:, index] *= values[i, 0]
                mesh._vpv[:, index] *= values[i, 1]
                mesh._vph[:, index] *= values[i, 2]
                mesh._vsv[:, index] *= values[i, 3]
                mesh._vsh[:, index] *= values[i, 4]
                mesh._eta[:, index] *= values[i, 5]
                mesh._qmu[index] *= values[i, 6]
                mesh._qkappa[index] *= values[i, 7]
        return mesh

    def __add__(self, other):
        assert np.allclose(self._vrmin, other._vrmin)
        assert np.allclose(self._vrmax, other._vrmax)
        model = self.__copy__()
        model._rho += other._rho
        model._vpv += other._vpv
        model._vph += other._vph
        model._vsv += other._vsv
        model._vsh += other._vsh
        model._eta += other._eta
        model._qmu += other._qmu
        model._qkappa += other._qkappa
        return model

    def _set_all_layers(self, index, values, scalar_value=0.):
        """
        Args:
            index (int): index of layer to be set
            values (ndarray): values.shape = (4,)
        """
        assert values.shape == (4,)
        assert values.dtype == np.float64
        self._rho[:, index] = values
        self._vpv[:, index] = values
        self._vph[:, index] = values
        self._vsv[:, index] = values
        self._vsh[:, index] = values
        self._eta[:, index] = values
        self._qmu[index] = scalar_value
        self._qkappa[index] = scalar_value

    def _add_boundary(self, r: float):
        """Add a boundary at radius=r (km).

        Returns:
            SeismicModel with added boundary
        """
        model = self.__copy__()
        if r in self._vrmin or r in self._vrmax:
            return model
        index = bisect.bisect_right(self._vrmin, r)
        model._vrmin = np.insert(model._vrmin, index, r)
        model._vrmax = np.insert(model._vrmax, index-1, r)
        model._rho = np.insert(
            model._rho, index-1, model._rho[:,index-1], axis=1)
        model._vpv = np.insert(
            model._vpv, index-1, model._vpv[:,index-1], axis=1)
        model._vph = np.insert(
            model._vph, index-1, model._vph[:,index-1], axis=1)
        model._vsv = np.insert(
            model._vsv, index-1, model._vsv[:,index-1], axis=1)
        model._vsh = np.insert(
            model._vsh, index-1, model._vsh[:,index-1], axis=1)
        model._eta = np.insert(
            model._eta, index-1, model._eta[:,index-1], axis=1)
        model._qmu = np.insert(
            model._qmu, index-1, model._qmu[index-1])
        model._qkappa = np.insert(
            model._qkappa, index-1, model._qkappa[index-1])
        return model

    def __copy__(self):
        """Deep copy of SeismicModel."""
        return SeismicModel(
            np.array(self._vrmin, dtype=np.float64),
            np.array(self._vrmax, dtype=np.float64),
            np.array(self._rho, dtype=np.float64),
            np.array(self._vpv, dtype=np.float64),
            np.array(self._vph, dtype=np.float64),
            np.array(self._vsv, dtype=np.float64),
            np.array(self._vsh, dtype=np.float64),
            np.array(self._eta, dtype=np.float64),
            np.array(self._qmu, dtype=np.float64),
            np.array(self._qkappa, dtype=np.float64))
    
    def get_zone(self, r):
        return bisect.bisect_right(self._vrmin, r) - 1

    def evaluate(self, r, poly):
        x = r / 6371.
        return (poly[0] 
               + poly[1]*x
               + poly[2]*x**2
               + poly[3]*x**3)
    
    def get_values(self, dr=1):
        rs = np.linspace(0, 6371, int(6371/dr))
        values = {'rho':[self.evaluate(r, self._rho[:, self.get_zone(r)])
                  for r in rs],
                  'vpv':[self.evaluate(r, self._vpv[:, self.get_zone(r)])
                  for r in rs],
                  'vph':[self.evaluate(r, self._vph[:, self.get_zone(r)])
                  for r in rs],
                  'vsv':[self.evaluate(r, self._vsv[:, self.get_zone(r)])
                  for r in rs],
                  'vsh':[self.evaluate(r, self._vsh[:, self.get_zone(r)])
                  for r in rs],
                  'eta':[self.evaluate(r, self._eta[:, self.get_zone(r)])
                  for r in rs],
                  'qmu':[self._qmu[self.get_zone(r)]
                  for r in rs],
                  'qkappa':[self._qkappa[self.get_zone(r)]
                  for r in rs]}
        return rs, values
    
    def plot(self, ax=None, parameters=None, color='black'):
        rs, values = self.get_values(dr=1.)
        if ax == None:
            fig, ax = plt.subplots(1,1)
            ax.set_prop_cycle(None)
        else:
            fig = None
        if parameters is None:
            parameters = values.keys() - {'qmu', 'qkappa', 'eta'}
        for i,key in enumerate(parameters):
            ax.plot(values[key], rs, label=key, color=color)
        ax.set_ylim(0, 6371)
        ax.set(
            xlabel='Velocity (km/s)',
            ylabel='Radius (km)')
        ax.legend()
        return fig, ax

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        with open(path, 'rb') as f:
            model = pickle.load(f)
        return model

if __name__ == '__main__':
    prem = SeismicModel.prem()
    # model parameters
    types = [ParameterType.VSV, ParameterType.VSH]
    radii = np.array([3480., 3630.], dtype=np.float64)
    model_params = ModelParameters(types, radii)
    # mesh
    model, mesh = prem.boxcar_mesh(model_params)
    # multiply mesh with values
    values_dict = {
        ParameterType.VSV: -0.1,
        ParameterType.VSH: 0.1}
    values_mat = model_params.get_values_matrix(values_dict)
    mesh_ = mesh.multiply(model_params.get_nodes(), values_mat)
    model_ = model + mesh_
    # figure
    rs, values = model_.get_values()
    rs_prem, values_prem = prem.get_values()
    fig, ax = plt.subplots(1,1)
    ax.plot(values_prem, rs_prem, color='blue')
    ax.plot(values, rs, color='red')
    ax.set_ylim(0, 6371)
    ax.set(
        xlabel='Velocity (km/s)',
        ylabel='Radius (km)')
    plt.show()