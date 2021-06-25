from dsmpy._tish import parameters
from dsmpy.modelparameters import ModelParameters, ParameterType
import bisect
import numpy as np
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
        model_id (str): model identifier (e.g., 'prem')
    
    Attributes:

    """

    def __init__(
            self, vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, model_id,
            mesh_type=None, model_params=None, discontinuous=False):
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
        self._model_id = model_id
        self._mesh_type = mesh_type
        self._model_params = model_params
        self._discontinuous = discontinuous

    def __copy__(self):
        vrmin = np.array(self._vrmin)
        vrmax = np.array(self._vrmax)
        rho = np.array(self._rho)
        vpv = np.array(self._vpv)
        vph = np.array(self._vph)
        vsv = np.array(self._vsv)
        vsh = np.array(self._vsh)
        eta = np.array(self._eta)
        qmu = np.array(self._qmu)
        qkappa = np.array(self._qkappa)
        model_id = str(self._model_id)
        mesh_type = str(self._mesh_type)
        model_params = self._model_params
        discontinuous = self._discontinuous
        return self.__new__(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, model_id,
            mesh_type, model_params, discontinuous)

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
    def get_model_id(self):
        return self._model_id
    def get_model_params(self):
        return self._model_params

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
        ], dtype=np.float32)
        vrmax = np.array([
            1217.5, 3479.5, 3631, 5611, 5711,
            5961, 6161, 6251, 6336.6, 6351, 6371
        ], dtype=np.float32)
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
        vsv = np.array([
            [3.667865, -0.001345, -4.440915, 0],
            [0, 0, 0, 0],
            [8.018341, -1.349895, 0, 0],
            [12.213901, -18.573085, 24.557329, -12.728015],
            [20.208945, -15.895645, 0, 0],
            [17.71732, -13.50652, 0, 0],
            [15.207335, -11.053685, 0, 0], # TODO check 15.212335
            [5.7502, -1.2742, 0, 0],
            [5.970824, -1.499059, 0, 0],
            [3.85, 0, 0, 0],
            [3.46, 0, 0, 0]
        ], dtype=np.float32).T
        vsh = np.array([
            [3.667865, -0.001345, -4.440915, 0],
            [0, 0, 0, 0],
            [8.018341, -1.349895, 0, 0],
            [12.213901, -18.573085, 24.557329, -12.728015],
            [20.208945, -15.895645, 0, 0],
            [17.71732, -13.50652, 0, 0],
            [15.207335, -11.053685, 0, 0],
            [5.7502, -1.2742, 0, 0],
            [5.970824, -1.499059, 0, 0],
            [3.85, 0, 0, 0],
            [3.46, 0, 0, 0]
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
        qmu = np.array([
            84.6, -1, 312, 312, 312, 143,
            143, 80, 600, 600, 600
        ], dtype=np.float32)
        qkappa = np.array([
            1327.7, 57823, 57823, 57823, 57823, 57823,
            57823, 57823, 57823, 57823, 57823
        ], dtype=np.float32)
        model_id = 'ak135'
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, model_id)

    @classmethod
    def prem(cls):
        """Return the Preliminary Reference Earth Model (PREM).

        References:
            Dziewonski and Anderson (1981)
        """
        vrmin = np.array([
            0, 1221.5, 3480, 3630, 5600, 5701, 5771,
            5971, 6151, 6291, 6346.6, 6356
            ], dtype=np.float32)
        vrmax = np.array([
            1221.5, 3480, 3630, 5600, 5701, 5771,
            5971, 6151, 6291, 6346.6, 6356, 6371
            ], dtype=np.float32)
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
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
        ], dtype=np.float32).T
        qmu = np.array([
            84.6, 1e12, 312, 312, 312, 143,
            143, 143, 80, 600, 600, 600
        ], dtype=np.float32)
        qkappa = np.array([
            1327.7, 57823, 57823, 57823, 57823,
            57823, 57823, 57823, 57823, 57823, 57823, 57823
        ], dtype=np.float32)
        model_id = 'prem'
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, model_id)

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
        model_id = 'iasp91'
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa, model_id)

    @classmethod
    def ak135_prime(cls):
        """Return model AK135.

        References:
            Kennett et al. (1995)
        """
        model = cls.ak135()
        model = model._del_boundary(6251)

        model._rho[:, -3] = np.array([5.36, -1, 0, 0])
        return model

    def lininterp_mesh(self, model_parameters, discontinuous=False):
        """
        Args:
            nodes (ndarray): nodes of the boxcar mesh
        Returns:
            model (SeismicModel): copy of self with added nodes
            mesh (SeismicModel): mesh with boxcar polynomials
        """
        model = self.__copy__()

        model = model._add_boundary(model_parameters.get_nodes()[0])
        model = model._add_boundary(model_parameters.get_nodes()[-1])

        for r in (set(model._vrmin) - set(model_parameters.get_nodes())):
            if (r > model_parameters.get_nodes()[0]
                and r < model_parameters.get_nodes()[-1]):
                model = model._del_boundary(r)

        for node in (set(model_parameters.get_nodes()) - set(model._vrmin)):
            model = model._add_boundary(node)

        for node in model_parameters.get_nodes()[:-1]:
            izone = model.get_zone(node)
            r0 = model._vrmin[izone]
            r1 = model._vrmax[izone]
            for param_type in ParameterType.structure_types(): # model_parameters._types
                if param_type == ParameterType.RADIUS:
                    continue
                y0 = model.get_value_at(r0, param_type) # self
                if discontinuous:
                    y1 = model.get_value_at(r1-1e-5, param_type) # self
                else:
                    y1 = model.get_value_at(r1, param_type) # self

                x0 = r0 / 6371.
                x1 = r1 / 6371.
                lin_elem = self._lin_element(x0, x1, y0, y1)
                model.set_value(izone, param_type, lin_elem)

        model._nzone = model._rho.shape[1]
        model._mesh_type = 'lininterp'
        model._model_params = model_parameters
        model._discontinuous = discontinuous
        return model

    def boxcar_mesh(self, model_parameters):
        """Create a boxcar mesh.

        Args:
            nodes (ndarray): nodes of the boxcar mesh

        Returns:
            model (SeismicModel): copy of self with added nodes
            mesh (SeismicModel): mesh with boxcar polynomials

        """
        model = self.__copy__()
        for node in model_parameters.get_nodes():
            model = model._add_boundary(node)
        # for i in range(model._nzone):
        #     model._set_all_layers(i, np.zeros(4, dtype=np.float32), 0.)
        nodes = model_parameters.get_nodes()
        for i in range(model_parameters._n_nodes-1):
            indexes = (set(np.where(model._vrmin < nodes[i+1])[0])
                & set(np.where(model._vrmin >= nodes[i])[0]))
            for type in [
                ParameterType.VSH, ParameterType.VSV,
                ParameterType.VPH, ParameterType.VPV,
                ParameterType.RHO, ParameterType.QMU,
                ParameterType.QKAPPA, ParameterType.ETA
            ]:
                value_avg = model.compute_avg(indexes, type)
                value_at_top = model.get_value_at(
                    model._vrmax[max(indexes)], type)
                for izone in indexes:
                    if type not in [ParameterType.QMU, ParameterType.QKAPPA]:
                        model.set_value(
                            izone, type, np.array([value_at_top, 0, 0, 0]))
                    else:
                        model.set_value(
                            izone, type, value_at_top)
        model._nzone = model._rho.shape[1]
        model._mesh_type = 'boxcar'
        model._model_params = model_parameters
        return model

    def compute_avg(
            self, izones: list, type: ParameterType, n=2) -> float:
        """Return the average value for type in layers izones.
        """
        if len(izones) == 0:
            return 0
        avg = 0.
        for izone in izones:
            rmin = self._vrmin[izone]
            rmax = self._vrmax[izone]
            radii = np.linspace(rmin, rmax, n, endpoint=False)
            values = [self.get_value_at(r, type) for r in radii]
            avg += sum([self.get_value_at(r, type) for r in radii])
        return avg / (n * len(izones))

    def triangle_mesh(self, model_parameters):
        """Create a triangular mesh.

        Args:
            nodes (ndarray): nodes of the triangle mesh (defines the
            radii of the triangles' peaks).

        Returns:
            model (SeismicModel): copy of self with added nodes
            mesh (SeismicModel): mesh with boxcar polynomials
            
        """
        model = self.__copy__()
        nodes = []
        for r in model_parameters._nodes:
            model = model._add_boundary(r)
            nodes.append(r)

        mesh = model.__copy__()
        for i in range(mesh._nzone):
            mesh._set_all_layers(i, np.zeros(4, dtype=np.float32), 0.)
        # all but top layer
        for i in range(len(nodes)-2):
            indexes = (set(np.where(mesh._vrmin < nodes[i+1])[0])
                & set(np.where(mesh._vrmin >= nodes[i])[0]))
            a0 = -nodes[i] / (nodes[i+1] - nodes[i])
            a1 = 1 / (nodes[i+1] - nodes[i]) * 6371.
            for index in indexes:
                mesh._set_all_layers(
                    index, np.array([a0, a1, 0, 0], dtype=np.float32), 1.)
        # top layer
        indexes = (set(np.where(mesh._vrmin < nodes[-1])[0])
                & set(np.where(mesh._vrmin >= nodes[-2])[0]))
        a0 = nodes[-1] / (nodes[-1] - nodes[-2])
        a1 = -1 / (nodes[-1] - nodes[-2]) * 6371.
        for index in indexes:
            mesh._set_all_layers(
                index, np.array([a0, a1, 0, 0], dtype=np.float32), 1.)
        model._nzone = model._rho.shape[1]
        mesh._nzone = mesh._rho.shape[1]
        mesh._mesh_type = 'triangle'
        mesh._model_params = model_parameters
        model._model_params = model_parameters
        return model, mesh

    def multiply(self, values: np.ndarray):
        """Return a copy of self with the model perturbations in values
        added.

        Args:
            values (np.ndarray): perturbations to model parameters
                (see ModelParameter.get_values_matrix())

        Returns:
            SeismicModel: new mesh with added values

        Examples:
            >>> model_params = ModelParameters(
                    types=[ParameterType.VSH],
                    radii=[3480., 3680.],
                    mesh_type='boxcar')
            >>> model = SeismicModel.prem().boxcar_mesh(model_params)
            >>> values_dict = {ParameterType.VSH: [0.1]}
            >>> values = model_params.get_values_matrix(values_dict)
            >>> updated_model = model.multiply(values)

            >>> model_params = ModelParameters(
                    types=[ParameterType.VSH],
                    radii=[3480., 3680.],
                    mesh_type='boxcar')
            >>> model = SeismicModel.prem().boxcar_mesh(model_params)
            >>> values = np.random.rand(
                    model_params.get_n_grd_params(), 9)
            >>> updated_model = model.multiply(values)

        """
        if self._model_params is None:
            return None
        mesh = self.__copy__()
        for i in range(len(self._model_params.get_nodes())-1):
            indexes = (
                set(np.where(
                    self._vrmin < self._model_params.get_nodes()[i+1])[0])
                & set(np.where(
                    self._vrmin >= self._model_params.get_nodes()[i])[0]))
            # indexes = (set(np.where(mesh._vrmin < nodes[i+1])[0])
            #     & set(np.where(mesh._vrmin >= nodes[i])[0]))
            if self._mesh_type == 'boxcar':
                for index in indexes:
                    mesh._rho[0, index] += values[i, 0]
                    mesh._vpv[0, index] += values[i, 1]
                    mesh._vph[0, index] += values[i, 2]
                    mesh._vsv[0, index] += values[i, 3]
                    mesh._vsh[0, index] += values[i, 4]
                    mesh._eta[0, index] += values[i, 5]
                    mesh._qmu[index] += values[i, 6]
                    mesh._qkappa[index] += values[i, 7]
            elif self._mesh_type == 'triangle':
                if i == 0:
                    for index in indexes:
                        mesh._rho[:, index] *= values[i, 0]
                        mesh._vpv[:, index] *= values[i, 1]
                        mesh._vph[:, index] *= values[i, 2]
                        mesh._vsv[:, index] *= values[i, 3]
                        mesh._vsh[:, index] *= values[i, 4]
                        mesh._eta[:, index] *= values[i, 5]
                        mesh._qmu[index] *= values[i, 6]
                        mesh._qkappa[index] *= values[i, 7]
                elif i == len(nodes) - 2:
                    a_add = np.array([1., 0., 0., 0.])
                    a_mul = np.array([-1., -1., 0., 0.])
                    for index in indexes:
                        mesh._rho[:, index] *= values[i-1, 0]
                        mesh._vpv[:, index] *= values[i-1, 1]
                        mesh._vph[:, index] *= values[i-1, 2]
                        mesh._vsv[:, index] *= values[i-1, 3]
                        mesh._vsh[:, index] *= values[i-1, 4]
                        mesh._eta[:, index] *= values[i-1, 5]
                        mesh._qmu[index] *= values[i-1, 6]
                        mesh._qkappa[index] *= values[i-1, 7]
                else:
                    for index in indexes:
                        a_add = np.array([1., 0., 0., 0.])
                        a_mul = np.array([-1., -1., 0., 0.])
                        mesh._rho[:, index] = (
                            mesh._rho[:,index] * values[i, 0]
                            + (a_mul*mesh._rho[:,index] + a_add)
                            * values[i-1, 0])
                        mesh._vpv[:, index] = (
                            mesh._vpv[:,index] * values[i, 1]
                            + (a_mul*mesh._vpv[:,index] + a_add)
                            * values[i-1, 1])
                        mesh._vph[:, index] = (
                            mesh._vph[:,index] * values[i, 2]
                            + (a_mul*mesh._vph[:,index] + a_add)
                            * values[i-1, 2])
                        mesh._vsv[:, index] = (
                            mesh._vsv[:,index] * values[i, 3]
                            + (a_mul*mesh._vsv[:,index] + a_add)
                            * values[i-1, 3])
                        mesh._vsh[:, index] = (
                            mesh._vsh[:,index] * values[i, 4]
                            + (a_mul*mesh._vsh[:,index] + a_add)
                            * values[i-1, 4])
                        mesh._eta[:, index] = (
                            mesh._eta[:,index] * values[i, 5]
                            + (a_mul*mesh._eta[:,index] + a_add)
                            * values[i-1, 5])
                        mesh._qmu[index] *= values[i, 6]
                        mesh._qkappa[index] *= values[i, 7]
            elif self._mesh_type == 'lininterp':
                for index in indexes:
                    # TODO check if no bug
                    r0 = self._vrmin[index]
                    r1 = self._vrmax[index]
                    # r0 = mesh._vrmin[index]
                    # r1 = mesh._vrmax[index]
                    r0_p = r0 + values[2 * i + 1, 8]
                    r1_p = r1 + values[2 * i + 3, 8]

                    for p_type in ParameterType.structure_types():
                        itype = p_type.value

                        # y0 = self.get_value_at(
                        #     r0, p_type) + values[2*i+1, itype]
                        # y1 = self.get_value_at(
                        #     r1-1e-5, p_type) + values[2*i+2, itype]

                        # TODO check
                        y0 = self.get_value_at(
                            r0_p, p_type) + values[2 * i + 1, itype]
                        y1 = self.get_value_at(
                            r1_p - 1e-5, p_type) + values[2 * i + 2, itype]

                        x0 = r0_p / 6371.
                        x1 = r1_p / 6371.
                        mesh.set_value(
                            index, p_type,
                            self._lin_element(x0, x1, y0, y1))

                    mesh._vrmin[index] = r0_p
                    mesh._vrmin[index+1] = r1_p
                    mesh._vrmax[index-1] = r0_p
                    mesh._vrmax[index] = r1_p

                    # TODO make further tests
                    # if r1_p < r1:
                    #     mesh = mesh._add_boundary(r1)
                    #     izone = mesh.get_zone(r1_p)
                    #     for p_type in ParameterType.structure_types():
                    #         mesh.set_value(
                    #             izone, p_type,
                    #             self.get_value(izone-1, p_type))
            else:
                raise ValueError(
                    'Expect "boxcar", "triangle", or "lininterp"')
        return mesh

    def gradient_models(self):
        """Return a list of seismic models to compute the waveform
        gradients with respect to model parameters using central
        differences.

        Returns:
            list of SeismicModel: the first half of the list corresponds
            to positive perturbations, the second half corresponds
            to negative perturbations
            list of float: corresponding perturbation vector

        Examples:
            >>> from dsmpy.dsm import compute_models_parallel
            >>> model_params = ModelParameters(
            ...        types=[ParameterType.VSH],
            ...        radii=[3480., 3680.],
            ...        mesh_type='boxcar')
            >>> model = SeismicModel.prem().boxcar_mesh(model_params)
            >>> grad_models, dxs = model.gradient_models()
            >>> outputs = compute_models_parallel(
            ...     dataset, grad_models, tlen=1638.4,
            ...     nspc=256, sampling_hz=20, mode=0)
            >>> n_params = len(grad_models) // 2
            >>> n_evs = len(outputs[0])
            >>> waveform_grads = []
            >>> for i in range(n_params):
            ...     waveform_grad_i = []
            ...     for iev in range(len(outputs[0])):
            ...         outputs[i][iev].to_time_domain()
            ...         outputs[i + n_params][iev].to_time_domain()
            ...         waveform_grad_i_iev = (outputs[i][iev].us
            ...             - outputs[i + n_params][iev].us) / (dxs[i]
            ...             - dxs[i + n_params])
            ...         waveform_grad_i.append(waveform_grad_i_iev)
            ...     waveform_grads.append(waveform_grad_i)
            >>> _, itypes, igrds = model_params.get_free_all_indices()
            >>> types = [model_params.types[i] for i in itypes]
            >>> radii = [model_params.get_grd_params()[i]
            ...          for i in igrds]

        """
        EPS = 1e-2
        if self._model_params is None:
            return None
        indices, itypes, igrds = self._model_params.get_free_all_indices()
        models = []
        dxs = []
        for count, i in enumerate(indices + indices):
            values = np.zeros(self._model_params.get_shape_value_matrix())
            type = self._model_params.get_types()[itypes[i]]
            if type in {ParameterType.QMU, ParameterType.QKAPPA,
                        ParameterType.RADIUS}:
                dx = 100 * EPS
            else:
                dx = EPS
            if count >= len(indices):
                dx = -dx
            values[igrds[i], type] = dx
            model_updated = self.multiply(values)
            models.append(model_updated)
            dxs.append(dx)
        return models, dxs

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

    def _lin_element(self, r0, r1, y0, y1):
        a = (r1*y0 - r0*y1) / (r1 - r0)
        b = (y1 - y0) / (r1 - r0)
        return np.array([a, b, 0, 0])

    def _set_all_layers(self, index: int, values, scalar_value=0.):
        """
        Initialize all layers and all types to values and scalara_values.

        Args:
            index (int): index of layer to be set
            values (ndarray): values.shape = (4,)
            scalar_value (float): value for QMU and QKAPPA (default is 0)
        """
        assert values.shape == (4,)
        # assert values.dtype == np.float64
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

        # TODO make sure no problem
        model._nzone += 1

        return model
    
    def _del_boundary(self, r: float):
        """Delete the boundary at radius=r (km).

        Returns:
            SeismicModel with removed boundary

        """
        model = self.__copy__()
        if (r not in self._vrmin) and (r not in self._vrmax):
            return model
        index = bisect.bisect_left(self._vrmin, r)
        model._vrmin = np.delete(model._vrmin, index)
        model._vrmax = np.delete(model._vrmax, index-1)
        model._rho = np.delete(
            model._rho, index, axis=1)
        model._vpv = np.delete(
            model._vpv, index, axis=1)
        model._vph = np.delete(
            model._vph, index, axis=1)
        model._vsv = np.delete(
            model._vsv, index, axis=1)
        model._vsh = np.delete(
            model._vsh, index, axis=1)
        model._eta = np.delete(
            model._eta, index, axis=1)
        model._qmu = np.delete(
            model._qmu, index)
        model._qkappa = np.delete(
            model._qkappa, index)
        return model

    def __copy__(self):
        """Deep copy of SeismicModel."""
        return SeismicModel(
            np.array(self._vrmin, dtype=np.float32),
            np.array(self._vrmax, dtype=np.float32),
            np.array(self._rho, dtype=np.float32),
            np.array(self._vpv, dtype=np.float32),
            np.array(self._vph, dtype=np.float32),
            np.array(self._vsv, dtype=np.float32),
            np.array(self._vsh, dtype=np.float32),
            np.array(self._eta, dtype=np.float32),
            np.array(self._qmu, dtype=np.float32),
            np.array(self._qkappa, dtype=np.float32),
            self._model_id,
            self._mesh_type,
            self._model_params,
            self._discontinuous)
    
    def get_zone(self, r: float) -> int:
        """Return index of layer that contains radius r
        (left inclusive, right exclusive).
        """
        return bisect.bisect_right(self._vrmin, r) - 1

    def evaluate(self, r, poly):
        x = r / 6371.
        return (poly[0] 
               + poly[1]*x
               + poly[2]*x**2
               + poly[3]*x**3)
    
    def get_values(self, dr=1.) -> (np.ndarray, dict):
        """Return a dict with values for each ParameterType.

        Args:
            dr (float): radius increment in km (default: 1)

        Returns:
            ndarray: radii
            dict: values. Keys are of type ParameterType

        """
        rs = np.linspace(0, 6371, int(6371/dr))
        values = {ParameterType.RHO:
                    [self.evaluate(r, self._rho[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.VPV:
                    [self.evaluate(r, self._vpv[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.VPH:
                    [self.evaluate(r, self._vph[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.VSV:
                    [self.evaluate(r, self._vsv[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.VSH:
                    [self.evaluate(r, self._vsh[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.ETA:
                    [self.evaluate(r, self._eta[:, self.get_zone(r)])
                    for r in rs],
                  ParameterType.QMU:
                    [self._qmu[self.get_zone(r)]
                    for r in rs],
                  ParameterType.QKAPPA:
                    [self._qkappa[self.get_zone(r)]
                    for r in rs]}
        return rs, values

    def get_value_at(self, r: float, type: ParameterType) -> float:
        """Return value at radius r.

        Args:
            r (float): radius
            type(ParameterType): type (e.g., ParameterType.VSH)

        Returns:
            float: value for type at radius r

        """
        if type == ParameterType.RHO:
            v = self.evaluate(r, self._rho[:, self.get_zone(r)])
        elif type == ParameterType.VPV:
            v = self.evaluate(r, self._vpv[:, self.get_zone(r)])
        elif type == ParameterType.VPH:
            v = self.evaluate(r, self._vph[:, self.get_zone(r)])
        elif type == ParameterType.VSV:
            v = self.evaluate(r, self._vsv[:, self.get_zone(r)])
        elif type == ParameterType.VSH:
            v = self.evaluate(r, self._vsh[:, self.get_zone(r)])
        elif type == ParameterType.ETA:
            v = self.evaluate(r, self._eta[:, self.get_zone(r)])
        elif type == ParameterType.QMU:
            v = self._qmu[self.get_zone(r)]
        elif type == ParameterType.QKAPPA:
            v = self._qkappa[self.get_zone(r)]
        return v

    def set_value(self, izone: int, type: ParameterType, values: np.ndarray):
        """Set the polynomial coefficients for layer izone.

        Args:
            izone (int): index of the layer
            type (ParameterType): type (e.g., ParameterType.VSH)
            values (np.ndarray): 3-degree polynomial coefficient

        """
        if type == ParameterType.RHO:
            self._rho[:, izone] = values
        elif type == ParameterType.VPV:
            self._vpv[:, izone] = values
        elif type == ParameterType.VPH:
            self._vph[:, izone] = values
        elif type == ParameterType.VSV:
            self._vsv[:, izone] = values
        elif type == ParameterType.VSH:
            self._vsh[:, izone] = values
        elif type == ParameterType.ETA:
            self._eta[:, izone] = values
        elif type == ParameterType.QMU:
            self._qmu[izone] = float(values)
        elif type == ParameterType.QKAPPA:
            self._qkappa[izone] = float(values)

    def get_value(self, izone: int, type: ParameterType) -> np.ndarray:
        """Get the polynomial coefficients for layer izone.

        Args:
            izone (int): index of the layer
            type (ParameterType): type (e.g., ParameterType.VSH)

        Returns:
            np.ndarray: 3-degree polynomial coefficients

        """
        if type == ParameterType.RHO:
            return self._rho[:, izone]
        elif type == ParameterType.VPV:
            return self._vpv[:, izone]
        elif type == ParameterType.VPH:
            return self._vph[:, izone]
        elif type == ParameterType.VSV:
            return self._vsv[:, izone]
        elif type == ParameterType.VSH:
            return self._vsh[:, izone]
        elif type == ParameterType.ETA:
            return self._eta[:, izone]
        elif type == ParameterType.QMU:
            return self._qmu[izone]
        elif type == ParameterType.QKAPPA:
            return self._qkappa[izone]

    def get_perturbations_to(
            self, model_ref, types, in_percent=False,
            range_dict=None):
        perturbations = np.zeros(
            (self._model_params._n_grd_params*len(types)), dtype=np.float32)
        for igrd in range(self._model_params._n_grd_params):
            ri = self._model_params._nodes[igrd]
            for ipar, param_type in enumerate(types):
                if param_type == ParameterType.RADIUS:
                    izone = self.get_zone(ri)
                    v = model_ref._vrmin[izone]
                    dv = self._vrmin[izone] - v
                else:
                    v = model_ref.get_value_at(ri, param_type)
                    dv = (self.get_value_at(ri, param_type)
                        - v)
                if in_percent:
                    dv /= v
                if range_dict is not None:
                    dv_min = range_dict[param_type][igrd, 0]
                    dv_max = range_dict[param_type][igrd, 1]
                    dv /= (dv_max - dv_min)
                index = igrd * len(types) + ipar
                perturbations[index] = dv
        # fig, ax = self.plot(types=[ParameterType.VSH])
        # model_ref.plot(types=[ParameterType.VSH], ax=ax)
        # plt.text(0, 4000, perturbations)
        # plt.show()
        return perturbations

    
    def plot(self, dr=1., ax=None, types=None, color=None, **kwargs):
        """Plot the seismicModel.

        Args:
            dr (float): depth increment (default is 1)
            ax (matplotlib.ax): ax
            types (ParameterTypes): e.g., RHO, VSH
            color (str): color
            kwargs (dict):

        Returns:
            Figure: matplotlib Figure object
            Axes: matplotlib Axes object

        """
        rs, values = self.get_values(dr=dr)
        if ax == None:
            fig, ax = plt.subplots(1,1)
            ax.set_prop_cycle(None)
        else:
            fig = None
        if types is None:
            types = (values.keys()
                - {ParameterType.QMU, ParameterType.QKAPPA, ParameterType.ETA})
        types = set(types) - {ParameterType.RADIUS}
        if 'label' in kwargs:
            label_kw = kwargs['label'] + ' '
            kwargs.pop('label', None)
        else:
            label_kw = ''
        for i,key in enumerate(types):
            label = label_kw + key.name
            if color is None:
                ax.plot(values[key], rs, label=label, **kwargs)
            else:
                ax.plot(values[key], rs, label=label, color=color, **kwargs)
        ax.set_ylim(0, 6371)
        ax.set(
            xlabel='Velocity (km/s)',
            ylabel='Radius (km)')
        # ax.legend()
        return fig, ax

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    
    def build_model(
            self, model, model_params,
            value_dict: dict):
        """Convenience function to build an updated seismic model with model
        perturbations added.

        Args:
            model (SeismicModel): the reference seismic model
            model_params (ModelParameters): model parameters
            value_dict_p (dict): dict of ParameterType:ndarray

        Returns:
            SeismicModel: the updated seismic model

        """
        values_mat = model_params.get_values_matrix(value_dict)
        model_updated = model.multiply(values_mat)
        return model_updated

    @staticmethod
    def load(path):
        with open(path, 'rb') as f:
            model = pickle.load(f)
        return model

if __name__ == '__main__':
    ak135 = SeismicModel.ak135()

    # model parameters
    types = [ParameterType.VSH]
    depth_moho = 6371. - 6336.6
    depth_410 = 410.
    depth_660 = 660.
    depth_max = 900. #1000
    
    n_upper_mantle = 0 #20
    n_mtz = 5 #10
    n_lower_mantle = 2 #12

    rs_upper_mantle = np.linspace(depth_410, depth_moho, n_upper_mantle)
    rs_mtz = np.linspace(depth_660, depth_410, n_mtz,
        endpoint=(n_upper_mantle==0))
    rs_lower_mantle = np.linspace(
        depth_max, depth_660, n_lower_mantle,
        endpoint=(n_mtz==0))
    radii = 6371. - np.round(
        np.hstack((rs_lower_mantle, rs_mtz, rs_upper_mantle)), 4)
    # print('dr_um={}, dr_mtz={}, dr_lm={}'.format(
    #     rs_upper_mantle[1] - rs_upper_mantle[0],
    #     rs_mtz[1] - rs_mtz[0],
    #     rs_lower_mantle[1] - rs_lower_mantle[0]))

    # radii = np.insert(radii, 2, 5700.)
    model_params = ModelParameters(types, radii, mesh_type='lininterp')

    model_ = ak135.lininterp_mesh(model_params, discontinuous=True)
    
    values_m = np.array(
        [0. if i%2==1 else 0. for i in range(model_params._n_grd_params)])
    values = np.array(
        [0. if i%2==1 else 0. for i in range(model_params._n_grd_params)])
    values[2] = -0.1
    values[3] = 0
    values_r = np.zeros(model_params._n_grd_params, dtype=np.float32)
    values_r[2] = -20
    values_r[3] = -0
    values_dict = {
        ParameterType.VSH: values,
        ParameterType.RADIUS: values_r}
    values_dict_m = {
        ParameterType.VSH: values_m,
        ParameterType.RADIUS: values_r}
    values_mat = model_params.get_values_matrix(values_dict)
    values_mat_m = model_params.get_values_matrix(values_dict_m)
    model_ = model_.multiply(model_params.get_nodes(), values_mat, values_mat_m)

    # # mesh
    # model, mesh = ak135.boxcar_mesh(model_params)

    # # multiply mesh with values
    # values = np.array(
    #     [0.1 * (-1)**i for i in range(model_params._n_grd_params)])
    # values_dict = {
    #     ParameterType.VSH: values}
    # values_mat = model_params.get_values_matrix(values_dict)
    # mesh_ = mesh.multiply(model_params.get_nodes(), values_mat)
    # model_ = model + mesh_

    # # retrieve model perturbations
    # perturbations = model_.get_perturbations_to(model, types)
    # print('Perturbations:', perturbations)

    # figure
    fig, ax = model_.plot(types=[ParameterType.VSH])
    ak135.plot(types=[ParameterType.VSH], ax=ax)
    ax.set_ylim([radii[0]-200, 6371.])
    plt.show()