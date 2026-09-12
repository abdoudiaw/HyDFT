"""Stanton & Murillo transport coefficients [PRE 93, 043203 (2016)] and the electron
viscosity of the Sci. Rep. 7, 15352 (2017) supplement.  Python mirror of
src/transport/stanton_murillo.f90 (same equations, same coefficients).

Reduced units (a = ion-sphere radius, wp = plasma frequency, k_B = 1):
    diffusion(Gamma, kappa)     D*   = D/(wp a^2)              Eq. 56
    viscosity(Gamma, kappa)     eta* = eta/(m n wp a^2)        Eq. 75
    conductivity(Gamma, kappa)  K*   = K/(n wp a^2)            Eq. 82
    mixture(Gamma12, kappa, x1, Z1, Z2, m1, m2) -> D12*, eta_tot*, K_tot*, k_T   Eqs. B8-B16
    electron_viscosity(rs, Gamma) -> eta_l*, eta_CV*, eta_SM*  supplement Eqs. 6-8

Physical (cgs) wrappers with the interface of the old SM.py, n [1/cm^3], m [g], Z, T [eV]:
    D_cgs, eta_cgs, K_cgs (single species), Dij_cgs, eta_tot_cgs, K_tot_cgs (binary).

Corrections with respect to the old SM.py:
  * K_13 b_1 = -0.38459 (was -0.384596) and K_22 a_1 = 0.85401 (was 0.8541), Table IV;
  * the single-species ionic screening term 4 pi Z^2 e^2 n/T had no Z^2 (lam_eff1);
  * with kappa = -1 the ion screening was counted twice (once inside lam_eff1, once
    again in Gamma sqrt(kappa^2 + 3 Gamma/(1 + 3 Gamma)));
  * like-species collision integrals Omega_ii used sqrt(2 pi/m_i) (eta2) or the
    inter-species reduced mass (Ktherm2) instead of mu_ii = m_i/2, i.e. sqrt(4 pi/m_i);
  * duplicated/undefined helpers (mfp, Lmf, sm.D) removed.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

SQRT3PI = np.sqrt(3 * np.pi)

# PRE 93, 043203 Table IV: (n,m) = (1,1), (1,2), (1,3), (2,2)
_A = {(1, 1): [1.4660, -1.7836, 1.4313, -0.55833, 0.061162],
      (1, 2): [0.52094, 0.25153, -1.1337, 1.2155, -0.43784],
      (1, 3): [0.30346, 0.23739, -0.62167, 0.56110, -0.18046],
      (2, 2): [0.85401, -0.22898, -0.60059, 0.80591, -0.30555]}
_B = {(1, 1): [0.081033, -0.091336, 0.051760, -0.50026, 0.17044],
      (1, 2): [0.20572, -0.16536, 0.061572, -0.12770, 0.066993],
      (1, 3): [0.68375, -0.38459, 0.10711, 0.10649, 0.028760],
      (2, 2): [0.43475, -0.21147, 0.11116, 0.19665, 0.15195]}


def Knm(n, m, g):
    """Reduced collision integral K_nm(g), fits (C22)-(C24)."""
    g = np.asarray(g, dtype=float)
    a, b = _A[(n, m)], _B[(n, m)]
    fact = float(math.factorial(m - 1))
    gw = np.minimum(np.maximum(g, 1e-300), 1.0)      # the weak-coupling polynomial is only used for g < 1
    weak = -0.25 * n * fact * np.log(sum(ak * gw**(k + 1) for k, ak in enumerate(a)))
    lg = np.log(np.maximum(g, 1e-300))
    strong = (b[0] + b[1] * lg + b[2] * lg**2) / (1 + b[3] * g + b[4] * g**2)
    return np.where(g < 1.0, weak, strong)


def K11(g): return Knm(1, 1, g)
def K12(g): return Knm(1, 2, g)
def K13(g): return Knm(1, 3, g)
def K22(g): return Knm(2, 2, g)


def g_eff(Gamma, kappa):
    """Plasma parameter with the effective screening length (28): Eq. 54."""
    return Gamma * np.sqrt(kappa**2 + 3 * Gamma / (1 + 3 * Gamma))


def g_mixture_factor(kappa, x, z, gkk):
    """a_tot/lambda_eff for a mixture (Eq. 68 with x_k in the numerator, see the Fortran header):
    x = number fractions, z = charge fractions Z_k n_k / sum Z_j n_j, gkk = Z_k^2 e^2/(a_tot T)."""
    s = kappa**2
    for xk, zk, gk in zip(x, z, gkk):
        if xk > 0:
            s += 3 * xk * gk / (1 + 3 * (xk / zk)**(1 / 3) * gk)
    return np.sqrt(s)


def diffusion(Gamma, kappa):
    return SQRT3PI / (12 * Gamma**2.5 * K11(g_eff(Gamma, kappa)))


def viscosity(Gamma, kappa):
    return 5 * SQRT3PI / (36 * Gamma**2.5 * K22(g_eff(Gamma, kappa)))


def conductivity(Gamma, kappa):
    return 25 * SQRT3PI / (48 * Gamma**2.5 * K22(g_eff(Gamma, kappa)))


def mixture(Gamma12, kappa, x1, Z1, Z2, m1, m2):
    """Binary mixture, first-order Chapman-Enskog (B8)-(B16), in units a_tot = 1 and the
    hydrodynamic plasma frequency w_hp = sqrt(4 pi <Z>^2 e^2 n/<m>):
    returns D12/(w_hp a^2), eta_tot/(rho w_hp a^2), K_tot/(n w_hp a^2), k_T."""
    x2 = 1.0 - x1
    T = Z1 * Z2 / Gamma12                     # e^2 = 1, a_tot = 1
    n = 3 / (4 * np.pi)
    rhoc = x1 * Z1 + x2 * Z2
    s = g_mixture_factor(kappa, [x1, x2], [x1 * Z1 / rhoc, x2 * Z2 / rhoc], [Z1**2 / T, Z2**2 / T])
    g11, g22, g12 = Z1**2 / T * s, Z2**2 / T * s, Z1 * Z2 / T * s
    mu12 = m1 * m2 / (m1 + m2)
    M1, M2 = m1 / (m1 + m2), m2 / (m1 + m2)
    om = lambda mu, zz, k: np.sqrt(2 * np.pi / (mu * T**3)) * zz**2 * k
    O11_22 = om(m1 / 2, Z1 * Z1, K22(g11))
    O22_22 = om(m2 / 2, Z2 * Z2, K22(g22))
    O12_11, O12_12, O12_13, O12_22 = (om(mu12, Z1 * Z2, K(g12)) for K in (K11, K12, K13, K22))
    eta1, eta2 = 5 * T / (8 * O11_22), 5 * T / (8 * O22_22)
    k1, k2 = 75 * T / (32 * m1 * O11_22), 75 * T / (32 * m2 * O22_22)
    A = O12_22 / (5 * O12_11)
    B = (5 * O12_12 - O12_13) / (5 * O12_11)
    C = 2 * O12_12 / (5 * O12_11) - 1
    E = T / (8 * M1 * M2 * O12_11)
    P1, P2 = O11_22 / (5 * M2 * O12_11), O22_22 / (5 * M1 * O12_11)
    Q1 = P1 * (6 * M2**2 + 5 * M1**2 - 4 * M1**2 * B + 8 * M1 * M2 * A)
    Q2 = P2 * (6 * M1**2 + 5 * M2**2 - 4 * M2**2 * B + 8 * M1 * M2 * A)
    Q12p = 15 * E * (P1 + P2 + (11 - 4 * B - 8 * A) * M1 * M2) / (2 * (m1 + m2))
    Q12 = 2 * P1 * P2 + 3 * (M1 - M2)**2 * (5 - 4 * B) + 4 * M1 * M2 * A * (11 - 4 * B)
    R1, R2 = 2 / 3 + M1 / M2 * A, 2 / 3 + M2 / M1 * A
    R12 = 4 * A / (3 * M1 * M2 * E) + E / (2 * eta1 * eta2)
    R12p = 4 / 3 + E / (2 * eta1) + E / (2 * eta2) - 2 * A
    S1 = M1 * P1 - M2 * (3 * (M2 - M1) + 4 * M1 * A)
    S2 = M2 * P2 - M1 * (3 * (M1 - M2) + 4 * M2 * A)
    den = x1**2 * Q1 + x2**2 * Q2 + x1 * x2 * Q12
    D12 = 3 * T / (16 * n * mu12 * O12_11)
    eta = (x1**2 * R1 + x2**2 * R2 + x1 * x2 * R12p) / (x1**2 * R1 / eta1 + x2**2 * R2 / eta2 + x1 * x2 * R12)
    K = (x1**2 * Q1 * k1 + x2**2 * Q2 * k2 + x1 * x2 * Q12p) / den
    kT = 5 * x1 * x2 * C * (x1 * S1 - x2 * S2) / den
    zbar, mbar = x1 * Z1 + x2 * Z2, x1 * m1 + x2 * m2
    whp = np.sqrt(4 * np.pi * zbar**2 * n / mbar)
    return D12 / whp, eta / (n * mbar * whp), K / (n * whp), kT


# ---------------------------------------------------------------- electrons (Sci. Rep. supplement)

def fermi_dirac(p, alpha):
    """I_p(alpha) = int_0^inf x^p/(1 + exp(x - alpha)) dx."""
    return quad(lambda x: x**p / (1.0 + np.exp(x - alpha)), 0.0, max(alpha, 0.0) + 80.0, limit=400)[0]


def electron_units(rs, Gamma):
    """theta = T/T_F and lambda_TF/a = 1/sqrt(4 pi e^2 dn/dmu) in HyDFT's reduced units
    (a = 1, m = 1, wp = 1: n0 = 3/4pi, kT = 1/(3 Gamma), e^2 = 1/3, hbar = 1/sqrt(3 rs))."""
    n0, beta, hbar, e2 = 3 / (4 * np.pi), 3.0 * Gamma, 1 / np.sqrt(3.0 * rs), 1 / 3.0
    cn = np.sqrt(2.0) / (np.pi**2 * hbar**3 * beta**1.5)
    alpha0 = brentq(lambda a: cn * fermi_dirac(0.5, a) - n0, -60.0, 400.0)
    tf = hbar**2 * (3 * np.pi**2 * n0)**(2 / 3) / 2.0
    theta = 1.0 / (beta * tf)
    dn_dmu = beta * cn * 0.5 * fermi_dirac(-0.5, alpha0)
    lam_tf = 1.0 / np.sqrt(4 * np.pi * e2 * dn_dmu)
    return dict(alpha0=alpha0, theta=theta, lambda_tf=lam_tf, hbar=hbar, beta=beta)


def eta_conti_vignale(rs):
    return 1 / np.sqrt(3 * rs) / (60 * rs**-1.5 + 80 / rs - 40 * rs**(-2 / 3) + 62 * rs**(-1 / 3))


def electron_viscosity(rs, Gamma):
    """eta_l* = (4/3)(eta_CV + theta eta_SM)/(1 + theta), eta_SM at g = Gamma a/lambda_TF.
    Returns (eta_l, eta_CV, eta_SM); the code's shear input is eta = 3 eta_l/4 (xi = 0)."""
    u = electron_units(rs, Gamma)
    sm = 5 * SQRT3PI / (36 * Gamma**2.5 * K22(Gamma / u['lambda_tf']))
    cv = eta_conti_vignale(rs)
    return 4 / 3 * (cv + u['theta'] * sm) / (1 + u['theta']), cv, sm


# ---------------------------------------------------------------- physical (cgs) wrappers

E2 = 1.44e-7          # e^2 in eV cm
ERG_PER_EV = 1.0 / 6.2415e11
HBAR = 6.5821e-16     # eV s
M_E = 9.109e-28       # g
KB = 1.380649e-16     # erg/K


def _wigner_seitz(n):
    return (3.0 / (4.0 * np.pi * n))**(1.0 / 3.0)


def _wp(n, Z, m):
    return np.sqrt(4 * np.pi * Z**2 * n * E2 * ERG_PER_EV / m)          # 1/s


def lambda_e(ne, T):
    """Electron Thomas-Fermi screening length [cm], Eq. 25: 1/lambda^2 = 4 pi e^2 n_e/sqrt(T^2 + (2 E_F/3)^2)."""
    EF = HBAR**2 * (3 * np.pi**2 * ne)**(2 / 3) / (2 * M_E) * ERG_PER_EV      # eV (hbar^2 k^2/m_e is in eV^2/erg)
    return 1.0 / np.sqrt(4 * np.pi * E2 * ne / np.sqrt(T**2 + (2 * EF / 3)**2))


def _single(n, m, Z, T, kappa):
    a = _wigner_seitz(n)
    Gamma = Z**2 * E2 / (a * T)
    if kappa is None:
        kappa = a / lambda_e(Z * n, T)
    return a, Gamma, kappa


def D_cgs(n, m, Z, T, kappa=None):
    """Self-diffusion [cm^2/s] and D*; kappa = a/lambda_e defaults to the TF electron screening."""
    a, Gamma, kappa = _single(n, m, Z, T, kappa)
    Ds = diffusion(Gamma, kappa)
    return Ds * _wp(n, Z, m) * a**2, Ds


def eta_cgs(n, m, Z, T, kappa=None):
    """Shear viscosity [g/(cm s)] and eta*."""
    a, Gamma, kappa = _single(n, m, Z, T, kappa)
    es = viscosity(Gamma, kappa)
    return es * m * n * _wp(n, Z, m) * a**2, es


def K_cgs(n, m, Z, T, kappa=None):
    """Thermal conductivity [erg/(cm s K)] (k_B included) and K* = K/(n k_B wp a^2)."""
    a, Gamma, kappa = _single(n, m, Z, T, kappa)
    ks = conductivity(Gamma, kappa)
    return ks * n * KB * _wp(n, Z, m) * a**2, ks


def _binary(n1, n2, m1, m2, Z1, Z2, T, kappa):
    n = n1 + n2
    a = _wigner_seitz(n)
    Gamma12 = Z1 * Z2 * E2 / (a * T)
    if kappa is None:
        kappa = a / lambda_e(Z1 * n1 + Z2 * n2, T)
    x1 = n1 / n
    zbar, mbar = x1 * Z1 + (1 - x1) * Z2, x1 * m1 + (1 - x1) * m2
    whp = np.sqrt(4 * np.pi * zbar**2 * n * E2 * ERG_PER_EV / mbar)
    return n, a, Gamma12, kappa, x1, mbar, whp


def Dij_cgs(n1, n2, m1, m2, Z1, Z2, T, kappa=None):
    """Interdiffusion coefficient [cm^2/s] of a binary mixture."""
    n, a, G, kappa, x1, mbar, whp = _binary(n1, n2, m1, m2, Z1, Z2, T, kappa)
    return mixture(G, kappa, x1, Z1, Z2, m1, m2)[0] * whp * a**2


def eta_tot_cgs(n1, n2, m1, m2, Z1, Z2, T, kappa=None):
    """Total viscosity [g/(cm s)] of a binary mixture."""
    n, a, G, kappa, x1, mbar, whp = _binary(n1, n2, m1, m2, Z1, Z2, T, kappa)
    return mixture(G, kappa, x1, Z1, Z2, m1, m2)[1] * n * mbar * whp * a**2


def K_tot_cgs(n1, n2, m1, m2, Z1, Z2, T, kappa=None):
    """Total thermal conductivity [erg/(cm s K)] of a binary mixture."""
    n, a, G, kappa, x1, mbar, whp = _binary(n1, n2, m1, m2, Z1, Z2, T, kappa)
    return mixture(G, kappa, x1, Z1, Z2, m1, m2)[2] * n * KB * whp * a**2


if __name__ == '__main__':
    print('K22(1) =', K22(1.0), ' D*(1, 1) =', diffusion(1.0, 1.0), '(MD: 2.54)')
    print('eta_l(rs=1.86, Gamma=1) =', electron_viscosity(1.86, 1.0))
    print('H-He 50/50, Gamma12 = 0.8, kappa = 0 (Fortran test prints the same case):', mixture(0.8, 0.0, 0.5, 1.0, 2.0, 1.0, 4.0))
