"""Generative and coordination functions of the vowel-space model.

Implements the three vocal-tract models of the article:

* the **generic model** (Sec. II-III): area functions made of the two
  first odd cosine components v1 = cos(pi x), v2 = cos(3 pi x) with
  coefficients (a1(theta), a2(theta)) given by the three-phase mixing
  function, then soft-rectified (`expif`, Eq. 4b);
* the **DRM** 4-tube model (Sec. V, Table II) driven by the coordination
  function with the {i,j,k} -> {Omega, Psi1, Psi2} transformation;
* **Fant's model** 4-parameter version (Table III): constriction position
  Xc, constriction area Ac, lip area Al and length L, embedded in the
  coordination function.

Everything follows the MATLAB supplementary scripts
(``simuGEN.m``, ``simuDRM.m``, ``simuFANT.m``, F. Berthommier 11/24/2021).
"""

import numpy as np


# ---------------------------------------------------------------- soft rect
def expif(vecin):
    """Soft rectifier, Eq. 4b (named ``expif``/``expcos`` in the MATLAB code).

    expif(y) = exp(y - 1) if y < 1 else y   -- C1 continuity at y = 1.
    """
    return np.where(vecin < 1, np.exp(vecin - 1), vecin)


#coordination transform
def ijk2psi(ijk):
    """{i, j, k} -> {Omega, Psi1, Psi2} transformation (Eq. 6 + selection).

    ``ijk`` rows are parameters, columns the {i, j, k} values attached to
    the phases theta_h in {pi/3, pi, 5 pi/3}.  For each row:
    Psi2 = atan(((i-k) sin(pi/3)) / (1/2 (i+k) - j)), then Psi1 is obtained
    from the first component h with h != Omega:
    Psi1 = (h - Omega) / cos(Psi2 - theta_h).
    """
    ijk = np.atleast_2d(np.asarray(ijk, dtype=float))
    a, b = ijk.shape
    om = np.zeros(a)
    psi1 = np.zeros(a)
    psi2 = np.zeros(a)
    phase = np.array([np.pi / 3, np.pi, 5 * np.pi / 3])
    for k in range(a):
        om[k] = np.mean(ijk[k, :])
        P1 = (ijk[k, 0] - ijk[k, 2]) * np.sin(np.pi / 3)
        P2 = 0.5 * (ijk[k, 0] + ijk[k, 2]) - ijk[k, 1]
        with np.errstate(divide='ignore', invalid='ignore'):
            psi2[k] = np.arctan(P1 / P2)            # atan(+/-inf) = +/- pi/2
        ind = np.where(ijk[k, :] - om[k])[0]        # first non-zero element
        psi1[k] = (ijk[k, ind[0]] - om[k]) / np.cos(psi2[k] - phase[ind[0]])
    return om, psi1, psi2


def mround(x):
    """MATLAB round: half away from zero (np.round is half-to-even)."""
    return np.floor(np.asarray(x) + 0.5).astype(int)


# ---------------------------------------------------------------- generic
ALPHA = 4.0 / 3.0 * np.sin(np.pi / 3)      # a2 amplitude of [o] and [e]

# characteristic phases and (a1, a2) of the 8 vowels [1,u,o,O,a,E,e,i]
VOWELS = ['ɨ', 'u', 'o', 'ɔ', 'a', 'ɛ', 'e', 'i']
INDTHETAB = np.array([1, 21, 31, 41, 61, 81, 91, 101])   # 1-based, step pi/60


def generic_P(theta, rho=1.0, n=100):
    """Three-phase mixing function P(x, theta), Eq. 3b (before rectification).

    P = 1 + 2 rho cos(theta) v1 + (4/3) rho sin(pi/3) sin(theta) v2,
    sampled at the n tubelet centers x = (i - 1/2)/n.
    """
    x = (np.arange(n) + 0.5) / n
    v1 = np.cos(np.pi * x)
    v2 = np.cos(3 * np.pi * x)
    return 1 + 2 * rho * np.cos(theta) * v1 + (4.0 / 3.0) * rho * np.sin(np.pi / 3) * np.sin(theta) * v2


def generic_area(P, L=17.5, n=None):
    """Rectified area function -> column format [length, area] for TLM3spec."""
    if n is None:
        n = len(P)
    return np.column_stack([np.full(n, L / n), expif(P)])


# ---------------------------------------------------------------- DRM
def drm_setup():
    """Coordination parameters {Omega, Psi1, Psi2} of the DRM (Table II).

    The {i, j, k} values are taken from the generic model at the cut
    points x in {0, L/3} (P1, P2); P3, P4 follow by internal antisymmetry.
    """
    ijk = np.array([
        [3.0, -1.0, 1.0],
        [0.5, 0.0, 2.5],
        2 - np.array([0.5, 0.0, 2.5]),
        2 - np.array([3.0, -1.0, 1.0]),
    ])
    return ijk2psi(ijk)


def drm_P(Pval, n=120):
    """DRM 4-tube -> n-tubelet area function (n multiple of 6)."""
    return np.concatenate([Pval[0] * np.ones(n // 6),
                           Pval[1] * np.ones(n // 3),
                           Pval[2] * np.ones(n // 3),
                           Pval[3] * np.ones(n // 6)])


def drm_area(Pval, L=17.5, n=120):
    P = drm_P(Pval, n)
    return np.column_stack([np.full(n, L / n), expif(P)])


# ---------------------------------------------------------------- Fant
def fant_setup(n=200):
    """{Omega, Psi1, Psi2} of Fant's model (Table III).

    Geometrical parameters: A = 1 cm2, L = 17.5 cm, dL = 1.5 cm,
    Lc = 0.9 n, l = 0.3 n (in tubelets).  Returns (Omega, Psi1, Psi2, Lc, l, Lip, A).
    """
    Lc = n * 0.9
    Lip = n - Lc
    l = n * 0.3
    A = 1.0
    delta = Lc - l
    Psi1 = np.array([0.3 * delta, 2 * A, -A, 1.5])
    Psi2 = np.array([5 * np.pi / 3, np.pi, np.pi / 3, np.pi / 3])
    Omega = np.array([Lc / 2, -1.5 * A, A / 2, 17.5])
    return Omega, Psi1, Psi2, Lc, l, Lip, A


def fant_P_raw(Pval, Lc, l, Lip, A):
    """Fant 4-tube reconstruction BEFORE soft rectification (rounding/merging)."""
    L1 = mround(Pval[0] - (l / 2))
    L2 = int(Lc - l - L1)
    return np.concatenate([A * np.ones(L1), Pval[1] * np.ones(int(l)),
                           A * np.ones(L2), Pval[2] * np.ones(int(Lip))])


def fant_area(Pval, Lc, l, Lip, A, n=200):
    """Fant 4-tube reconstruction: rounding, merging, soft rectification."""
    P = fant_P_raw(Pval, Lc, l, Lip, A)
    return np.column_stack([np.full(n, Pval[3] / n), expif(P)])


# ---------------------------------------------------------------- Fourier
def fourier_a12(P, n=None):
    """Cosine Fourier coefficients (a1, a2) of an area function, Eq. 8.

    a1 = 2/n sum A(i) cos(pi (i-1/2)/n),  a2 = same with 3 pi.
    """
    P = np.asarray(P)
    if n is None:
        n = len(P)
    x = (np.arange(len(P)) + 0.5) / n
    v1 = np.cos(np.pi * x)
    v2 = np.cos(3 * np.pi * x)
    return 2 * np.sum(P * v1) / n, 2 * np.sum(P * v2) / n


def se_estimates(f1b, f2b, a1, a2):
    """Formant estimates from the coefficients (Schroeder-Ehrenfest + bias).

    Returns dict with:
      'est'  : df1 = -a1/2 - (a2/2)^2, df2 = -a2/2   (biased "est" relation)
      'se'   : df1 = -a1/2,           df2 = -a2/2   (SE from rectified coef.)
    """
    return {
        'est': (f1b - f1b * ((a1 / 2) + (a2 / 2) ** 2), f2b - f2b * a2 / 2),
        'se': (f1b - f1b * a1 / 2, f2b - f2b * a2 / 2),
    }
