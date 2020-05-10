

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

    @classmethod
    def prem(cls):
        """Return the Preliminary Reference Earth Model (PREM).

        References:
            Dziewonski and Anderson (1981)
        """
        vrmin = 
        vrmax =
        rho =
        vpv =
        vph = 
        vsv = 
        vsh = 
        eta = 
        qmu = 
        qkappa = 
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)

    @classmethod
    def ak135(cls):
        """Return model AK135.

        References:
            Kennett et al. (1995)
        """
        vrmin =
        vrmax =
        rho =
        vpv =
        vph = 
        vsv = 
        vsh = 
        eta = 
        qmu = 
        qkappa = 
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)

    @classmethod
    def iasp91(cls):
        """Return model IAS91.

        References:
        """
        vrmin =
        vrmax =
        rho =
        vpv =
        vph = 
        vsv = 
        vsh = 
        eta = 
        qmu = 
        qkappa = 
        return cls(
            vrmin, vrmax, rho, vpv, vph,
            vsv, vsh, eta, qmu, qkappa)