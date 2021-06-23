import numpy as np
from enum import IntEnum


class ModelParameters:
    """Represent parameters for a seismic model.

    The parameters are specified by:
        - their types: ParameterType.VSH, VSV, VPH, VPV, ETA, QMU,
            QKAPPA
        - their radii (in km): the model parameters represent a layered
            model, and the radii specify the boundaries of these layers.
        - the mesh type:
            - 'boxcar': constant properties in each layer
            - 'lininterp': discontinuous linear spline in each layer
            - 'triangle': still in development

    Args:
        types (list of ParameterType): e.g., RHO, VSH
        radii (ndarray): radii of boundaries of perturbed layers
        mesh_type (str): 'boxcar', 'triangle', or 'lininterp'

    """

    def __init__(self, types, radii, mesh_type='boxcar'):
        self._types = types
        self._radii = radii
        self._mesh_type = mesh_type
        self.mask_dict = None
        self.equal_dict = None
        self.discon_arr = None
        if mesh_type == 'boxcar':
            self._n_nodes = len(radii)
            self._n_grd_params = len(radii) - 1
            self._nodes = np.array(radii)
        elif mesh_type == 'triangle':
            # self._n_nodes = int(
            #     len(radii) - 1 + np.ceil((len(radii)-1) / 2))
            self._n_nodes = len(radii)
            self._n_grd_params = 2 * len(radii) - 1
            n = 2 * len(radii) - 1
            self._nodes = np.zeros(n, dtype=np.float32)
            for i in range(n - 1):
                self._nodes[i] = (radii[i // 2] + (i % 2)
                                  * (radii[i // 2 + 1] - radii[i // 2]) / 2.)
            self._nodes[-1] = radii[-1]
        elif mesh_type == 'lininterp':
            self._n_nodes = len(radii)
            self._n_grd_params = len(radii) * 2  # TODO check consequences
            self._nodes = np.array(radii)
        self.it = 0
        self.set_constraints()

    def get_nodes(self) -> np.ndarray:
        """Return an array of radial grid nodes."""
        return self._nodes

    def get_n_nodes(self) -> int:
        """Return the number of radial grid nodes."""
        return self._n_nodes

    def get_n_grd_params(self) -> int:
        """Get the number of radial parameters at which structure
        parameters (e.g. ParameterType.VSH) can be set.

        The number of parameters depend on the mesh type:
            - boxcar: get_n_grd_params() = get_nodes() - 1
            - lininterp: get_n_grd_params() = get_nodes() * 2.
              At each radial node, there is one parameter for the line
              segment above the node, and one parameter for the line
              segment below the node.
        """
        return self._n_grd_params

    def get_n_params(self) -> int:
        """Return the total number of parameters
        get_n_grd_params() * len(get_types())."""
        return self._n_grd_params * len(self._types)

    def get_mesh_type(self) -> int:
        """Return the mesh type."""
        return self._mesh_type

    def get_types(self) -> int:
        """Return the parameter types"""
        return self._types

    def set_constraints(
            self, mask_dict=None, equal_dict=None,
            discon_arr=None):
        """Set constraints for model parameters at each grid point.

        Three types of constraits. mask_dict and equal_dict act on
        the n_grd_params grid parameters; discon_arr acts on the n_nodes
        radial nodes.
            - mask_dict: fix the value for a given ParameterType
                and grid point (among n_grd_params)
            - equal_dict: for a given ParameterType, an integer array
                specify grid points that takes the same value
            - discon_arr: for each radial node, specify wether this node
                is discontinuous or not

        Args:
            mask_dict (dict): keys are of type ParameterType,
                values are boolean np.ndarray of shape (n_grd_params,)
            equal_dict (dict): keys are of type ParameterType,
                values are integer np.ndarray of shape (n_grd_params,)
            discon_arr (np.ndarray): boolean array of shape (n_nodes,)

        Examples:
            >>> mask_dict[ParameterType.VSH] = np.ones(
                    model_params.get_n_grd_params(), dtype='bool')
            >>> equal_dict[ParameterType.VSH] = np.ones(
                    model_params.get_n_grd_params(), dtype='bool')
            >>> discon_arr = np.ones(
                    model_params.get_n_nodes(), dtype='bool')
            >>> model_params.set_constraints(
                    mask_dict, equal_dict, discon_arr)

        """
        if mask_dict is not None:
            self.mask_dict = mask_dict
        else:
            self.mask_dict = dict()
            for param_type in self._types:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params, dtype='bool')
        if equal_dict is not None:
            self.equal_dict = equal_dict
        else:
            self.equal_dict = dict()
            for param_type in self._types:
                self.equal_dict[param_type] = np.arange(
                    self._n_grd_params, dtype='int')
        if discon_arr is not None:
            self.discon_arr = discon_arr
        else:
            self.discon_arr = np.ones(
                self._n_nodes, dtype='bool')

        for param_type in self._types:
            if param_type not in self.equal_dict:
                self.equal_dict[param_type] = np.arange(
                    self._n_grd_params, dtype='int')
            if param_type not in self.mask_dict:
                self.mask_dict[param_type] = np.ones(
                    self._n_grd_params, dtype='bool')

    def next(self) -> (int, int, int):
        """Increment and return the indices of the next model parameter.
        next() skips indices with constraints.

        Returns:
            int: iteration counter
            int: index of the current type (get_types()[i])
            int: index of the current radial parameter
        """
        found = False
        istop = 0
        seen_igrd = set()
        while not found and istop < 1000:
            self.it = self.it % self.get_n_params()
            itype = int(self.it // self._n_grd_params)
            igrd = (self.it % self._n_grd_params)
            p_type = self._types[itype]
            if (self.mask_dict[p_type][igrd]
                    and self.equal_dict[p_type][igrd] == igrd):
                found = True
            else:
                self.it += 1
            istop += 1

            if (igrd not in seen_igrd
                    and
                    (not self.discon_arr[igrd//2]
                     or p_type == ParameterType.RADIUS)):
                self.it += 1
                seen_igrd.add(igrd)
                igrd += 1
        return self.it, itype, igrd

    def get_free_indices(self) -> np.ndarray:
        """Return the indices of parameters without constraints.
        """
        indices_lst = []
        it_0 = self.it
        self.it = 0
        i = 0
        stop = False
        while not stop:
            i, _, _ = self.next()
            if i in indices_lst:
                stop = True
            else:
                indices_lst.append(i)
            self.it += 1
        self.it = it_0
        return np.array(indices_lst)

    def apply_constraint(self, values: np.ndarray) -> np.ndarray:
        """Apply the model parameter constraints to the valye matrix.

        Args:
            values (np.ndarray): (n_grd_params, 9)-matrix with model
                perturbations

        Returns:
            np.ndarray: a copy of values with constraints applied

        """
        values_mat = np.array(values)
        if self.mask_dict is not None:
            for key, mask in self.mask_dict.items():
                values_mat[~mask, key] = 0.
        if self.discon_arr is not None:
            for i in np.where(~self.discon_arr)[0]:
                indices_struct = ParameterType.structure_types()
                values_mat[2*i, indices_struct] = values_mat[
                    2*i+1, indices_struct]
        if self.equal_dict is not None:
            for key, indexes in self.equal_dict.items():
                for i, j in enumerate(indexes):
                    if i > j:
                        values_mat[i, key] = values_mat[j, key]
        return values_mat

    def get_values_matrix(
            self, values_dict: dict) -> np.ndarray:
        """Get the matrix used in Seismicmodel.multiply"""
        values_mat = np.zeros((self._n_grd_params, 9), np.float64)
        for key, values in values_dict.items():
            values_mat[:, key] = values

        # TODO check that it's done elsewhere, or modify code
        # TODO check that it doesn't change behavior
        if self.mask_dict is not None:
            for key, mask in self.mask_dict.items():
                values_mat[~mask, key] = 0.
        if self.discon_arr is not None:
            for i in np.where(~self.discon_arr)[0]:
                indices_struct = ParameterType.structure_types()
                values_mat[2*i, indices_struct] = values_mat[
                    2*i+1, indices_struct]
        if self.equal_dict is not None:
            for key, indexes in self.equal_dict.items():
                for i, j in enumerate(indexes):
                    if i > j:
                        values_mat[i, key] = values_mat[j, key]
        return values_mat


class ParameterType(IntEnum):
    RHO = 0
    VPV = 1
    VPH = 2
    VSV = 3
    VSH = 4
    ETA = 5
    QMU = 6
    QKAPPA = 7
    RADIUS = 8

    @staticmethod
    def structure_types() -> list:
        """Return a list of structural parameters:
        [RHO, VPV, VPH, VSV, VSH, ETA]"""
        return [ParameterType(i) for i in range(6)]


if __name__ == '__main__':
    types = [ParameterType.VSV, ParameterType.VSH]
    radii = np.array([3480., 3700.], dtype=np.float32)
    model_params = ModelParameters(types, radii)
    values_dict = {
        ParameterType.VSV: -1.,
        ParameterType.VSH: 1.}
    values_mat = model_params.get_values_matrix(values_dict)
    print(values_mat)
