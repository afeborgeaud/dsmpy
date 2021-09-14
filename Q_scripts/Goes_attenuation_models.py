import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from collections import namedtuple

"""
This script implements the mantle seismic attenuation model
of Goes et al. (2004) and Maguire et al. (2016).
Thanks to Saskia Goes for her detailed instructions
and for detailing her favoured model.

Bob Myhill, 2021/06/04
"""

AnelasticProperties = namedtuple('AnelasticProperties',
                                 ['V_P', 'V_S',
                                  'Q_S', 'Q_K', 'Q_P', 'T_solidus'])

"""
First, we load a couple of tables.
PREM is required for the pressure to depth conversion
The geotherm is required only for plotting purposes
"""
prem = np.loadtxt('prem.txt', unpack=True)
prem_depths = prem[0]
prem_pressures = prem[2]
prem_QK_L6 = prem[6]
prem_QS_L6 = prem[7]
fn_pressure = interp1d(prem_depths, prem_pressures, kind='linear')


# Here we fit lower mantle solidus curve to Fiquet et al., 2010,
# using a pin at the high pressure end of the Herzberg study.
# The fitted parameters can be inserted into the
# Simon_Glatzel_lower_mantle_Fiquet function, below.
if False:
    def Simon_Glatzel(P, a, b):
        Pr = 36.e9
        Tr = 2800.
        return Tr*np.power(((P - Pr)/a + 1.), b)

    Pm = np.array([22.5e9, 36.e9, 135.e9])
    Tm = np.array([1086. - 5.7*22.5 + 390*np.log(22.5) + 273.15,
                   2800.,
                   4180.])

    popt, pcov = curve_fit(Simon_Glatzel, Pm, Tm, [40.e9, 0.3])

def Simon_Glatzel_lower_mantle_Fiquet(P):
    """
    Simon Glatzel model to fit estimates of melting point
    at 36 and 135 GPa, pinned at the 22.5 GPa Herzberg estimate

    The values of a and b are calculated using the fitting
    procedure above (commented out)
    """
    Pr = 36.e9
    Tr = 2800.
    a = 3.86695953e+10
    b = 3.15554341e-01
    return Tr*np.power(((P - Pr)/a + 1.), b)

def peridotite_solidus(P):
    """
    Returns an estimate of the peridotite solidus, using three studies:
    Hirschmann (2000) (0 GPa, then linear extrapolation to 2.7 GPa)
    Herzberg et al (2000) (2.7 - 22.5 GPa)
    Fiquet et al. (2010) (> 22.5 GPa)

    This curve is continuous with pressure, but not differentiable.
    """
    if P < 2.7e9: # interpolation between Hirschmann (2000) at 0 GPa and Herzberg et al (2000) at 2.7GPa
        T = 1120.661 +  (1086. - 5.7*2.7 + 390*np.log(2.7) - 1120.661)*P/2.7e9 + 273.15
    elif P < 22.5e9: # Herzberg et al (2000)
        T = 1086. - 5.7*P/1.e9 + 390*np.log(P/1.e9) + 273.15
    else: # Fiquet et al. (2010)
        T = Simon_Glatzel_lower_mantle_Fiquet(P)
    return T

def anelastic_properties(elastic_Vp, elastic_Vs, depth, temperature,
                         frequency, Q_model, dT_Q_constant_above_solidus=0):
    """
    Calculates the anelastic Vp and Vs, QS, and QK according to the model
    used in Maguire et al., 2016.

    The effects of anelasticity on shear wave velocity are incorporated using
    a model for the S-wave quality factor QS that varies with
    temperature T and depth z as
    QS(z,T) = Qo ω a exp(a ξ Tm / T), where
    ω is frequency,
    a is exponential frequency dependence,
    ξ is a depth scaling factor and
    Tm is the dry solidus melting temperature.

    The solidus is provided by the function peridotite_solidus.
    The conversion from pressure to depth is from PREM (Dziewonski et al.)

    To avoid step-changes in QS at the top and base of
    the mantle, transition regions 2.2 GPa wide are implemented.
    At a reference temperature of 750K, the center of the ol-wd transition
    is at 11.1 GPa. At the same reference temperature, the center
    of the postspinel transition is at 26.1 GPa. Clapeyron slopes of
    2.4e6 Pa/K and -2.2e6 Pa/K are applied.

    QK is chosen to be temperature independent.

    The anelastic seismic velocities are then calculated as follows:
    lmda = 4/3 * (elastic_Vs/elastic_Vp)^2
    1/QP = (1. - lmda)/QK + lmda/QS

    If 1/QP is negative, it is set to 0.

    anelastic_Vp = elastic_Vp*(1 - invQP/(2tan(pi*alpha/2)))
    anelastic_Vs = elastic_Vs*(1 - invQS/(2tan(pi*alpha/2)))

    Arguments
    ---------
    elastic_Vp : the elastic P-wave velocity
    elastic_Vs : the elastic S-wave velocity
    depth : the depth in m
    temperature : the temperature in K
    frequency: the frequency of the seismic waves in Hz
    Q_model: a dictionary containing the parameters for the attenuation model
    dT_Q_constant_above_solidus: if the temperature > (solidus temperature + dT),
    the value of QS, QK and a are frozen at the values corresponding to
    (solidus temperature + dT).

    Returns
    -------

    An instance of an AnelasticProperties named tuple.
    Has the following attributes:
    V_P, V_S, Q_S, Q_K, Q_P
    """

    pressure = fn_pressure(depth)
    Tm = peridotite_solidus(pressure)

    # Freezes QS if above a certain temperature
    if dT_Q_constant_above_solidus < temperature - Tm:
        Q_temperature = Tm + dT_Q_constant_above_solidus
    else:
        Q_temperature = temperature

    QSs = []
    QKs = []
    alphas = []
    for Q_mod in [Q_model['UM'], Q_model['TZ'], Q_model['LM']]:
        QSs.append((Q_mod['Q0'] *
                    np.power(frequency, Q_mod['a']) *
                    np.exp(Q_mod['a']*Q_mod['g']*Tm/Q_temperature)))
        QKs.append(Q_mod['QK'])
        alphas.append(Q_mod['a'])

    P_smooth_halfwidth = 1.1e9
    T_ref = 750. # K
    pressure_tztop = 11.1e9 + 2.4e6*(temperature - T_ref)
    pressure_tzbase = 26.1e9 - 2.2e6*(temperature - T_ref)

    if pressure < pressure_tztop - P_smooth_halfwidth:
        QS = QSs[0]
        QK = QKs[0]
        alpha = alphas[0]

    elif pressure < pressure_tztop + P_smooth_halfwidth:
        f = (pressure - (pressure_tztop - P_smooth_halfwidth)) / (2.*P_smooth_halfwidth)
        QS = (1. - f)*QSs[0] + f*QSs[1]
        QK = (1. - f)*QKs[0] + f*QKs[1]
        alpha = (1. - f)*alphas[0] + f*alphas[1]

    elif pressure < pressure_tzbase - P_smooth_halfwidth:
        QS = QSs[1]
        QK = QKs[1]
        alpha = alphas[1]

    elif pressure < pressure_tzbase + P_smooth_halfwidth:
        f = (pressure - (pressure_tzbase - P_smooth_halfwidth)) / (2.*P_smooth_halfwidth)
        QS = (1. - f)*QSs[1] + f*QSs[2]
        QK = (1. - f)*QKs[1] + f*QKs[2]
        alpha = (1. - f)*alphas[1] + f*alphas[2]

    else:
        QS = QSs[2]
        QK = QKs[2]
        alpha = alphas[2]

    invQS = 1./QS
    invQK =1./QK

    lmda = 4./3.*np.power(elastic_Vs/elastic_Vp, 2.)
    invQP = (1. - lmda)*invQK + lmda*invQS

    if invQP < 0.:
        invQP = 0.
        QP = np.inf
    else:
        QP = 1./invQP

    anelastic_Vp = elastic_Vp*(1. - invQP/(2.*np.tan(np.pi*alpha/2.)))
    anelastic_Vs = elastic_Vs*(1. - invQS/(2.*np.tan(np.pi*alpha/2.)))

    return AnelasticProperties(V_P=anelastic_Vp,
                               V_S=anelastic_Vs,
                               Q_S=QS,
                               Q_K=QK,
                               Q_P=QP,
                               T_solidus=Tm)

# Q4g - low T dependence (after Goes et al. 2004)
Q4g = {'UM': {'Q0': 0.1,
              'g': 38.,
              'a': 0.15,
              'QK': 1000.},
       'TZ': {'Q0': 3.5,
              'g': 20.,
              'a': 0.15,
              'QK': 1000.},
       'LM': {'Q0': 35.,
              'g': 10.,
              'a': 0.15,
              'QK': 1000.}}

# Q6g - strong T dependence (after Goes et al. 2004
Q6g = {'UM': {'Q0': 0.1,
              'g': 38.,
              'a': 0.15,
              'QK': 1000.},
       'TZ': {'Q0': 0.5,
              'g': 30.,
              'a': 0.15,
              'QK': 1000.},
       'LM': {'Q0': 3.5,
              'g': 20.,
              'a': 0.15,
              'QK': 1000.}}

# Q7g - intermediate T dependence
# (most consistent with Matas and Bukuwinski 2007)
Q7g = {'UM': {'Q0': 0.1,
              'g': 38.,
              'a': 0.15,
              'QK': 1000.},
       'TZ': {'Q0': 0.5,
              'g': 30.,
              'a': 0.15,
              'QK': 1000.},
       'LM': {'Q0': 1.5,
              'g': 26.,
              'a': 0.15,
              'QK': 1000.}}


if __name__ == "__main__":

    geotherm = np.loadtxt('anderson_82.txt', unpack=True)
    depths = geotherm[0]
    temperatures = geotherm[1]
    fn_geotherm = interp1d(depths, temperatures, kind='linear')
    
    frequency = 1. # Hz

    depths = np.linspace(0., 2880.e3, 1001)
    Vps = np.empty_like(depths)
    Vss = np.empty_like(depths)
    QSs = np.empty_like(depths)
    QKs = np.empty_like(depths)
    Vpsh = np.empty_like(depths)
    Vssh = np.empty_like(depths)
    QSsh = np.empty_like(depths)
    QKsh = np.empty_like(depths)
    temperatures = np.empty_like(depths)
    melting_temperatures = np.empty_like(depths)
    pressures = np.empty_like(depths)

    dT = 300.
    for i, z in enumerate(depths):
        temperatures[i] = fn_geotherm(z)

        # modify the temperature profile to make it continuous with depth
        if z > 670.e3:
            temperatures[i] -= 150.

        pressures[i] = fn_pressure(z)
        melting_temperatures[i] = peridotite_solidus(pressures[i])

        p0 = anelastic_properties(elastic_Vp=1.,
                                  elastic_Vs=1.,
                                  depth=z,
                                  temperature=temperatures[i],
                                  frequency=frequency,
                                  Q_model=Q7g)

        Vps[i], Vss[i], QSs[i], QKs[i] = [p0.V_P, p0.V_S, p0.Q_S, p0.Q_K]

        p1 = anelastic_properties(elastic_Vp=1.,
                                  elastic_Vs=1.,
                                  depth=z,
                                  temperature=temperatures[i]+dT,
                                  frequency=frequency,
                                  Q_model=Q7g)

        Vpsh[i], Vssh[i], QSsh[i], QKsh[i] = [p1.V_P, p1.V_S, p1.Q_S, p1.Q_K]


    fig = plt.figure(figsize=(10,8))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

    ax[0].plot(depths/1.e3, Vps, label='Vp')
    ax[0].plot(depths/1.e3, Vss, label='Vs')
    ax[0].plot(depths/1.e3, Vpsh, label=f'Vp (+{dT} K)')
    ax[0].plot(depths/1.e3, Vssh, label=f'Vs (+{dT} K)')

    ax[1].plot(depths/1.e3, pressures/1.e9, label='PREM')

    #*; -150 K at z > 670 km
    ax[2].plot(depths/1.e3, temperatures, label='geotherm (Anderson, 1982*)')
    ax[2].plot(depths/1.e3, temperatures+300, label=f'geotherm (Anderson, 1982*) + {dT} K')
    ax[2].plot(depths/1.e3, melting_temperatures, label='solidus')

    ax[3].plot(depths/1.e3, QSs, label='QS')
    ax[3].plot(depths/1.e3, QKs, label='QK')

    ax[3].plot(depths/1.e3, QSsh, label=f'QS (+{dT} K)')
    ax[3].plot(depths/1.e3, QKsh, label=f'QK (+{dT} K)')

    ax[3].plot(prem_depths/1.e3, prem_QS_L6, label='QS (QL6; Durek and Ekstrom, 1996)')
    ax[3].plot(prem_depths/1.e3, prem_QK_L6, label='QK (QL6; Durek and Ekstrom, 1996)')

    ax[0].set_ylabel('Anelastic velocity / elastic velocity')
    ax[1].set_ylabel('Pressure (GPa)')
    ax[2].set_ylabel('Temperature (K)')
    ax[3].set_ylabel('Q')

    ax[3].set_ylim(1e1, 1.e6)
    ax[3].set_xlim(0, 2880.)
    for i in range(4):
        ax[i].set_xlabel('Depth (km)')
        ax[i].legend()

    ax[3].set_yscale('log')

    fig.tight_layout()

    fig.savefig('Q7g_model_predictions.pdf')
    plt.show()
