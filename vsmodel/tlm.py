"""Transmission-line model (TLM) of the vocal tract -- Python translation.

Faithful translation of the ICP-Grenoble MATLAB code by Pierre Badin
(``TLM.m`` / ``TLM3spec.m`` and their subfunctions, c. 1984-2006), as used in
the supplementary material of:

    F. Berthommier, "A mathematical model of the vowel space",
    JASA-Express Letters (2021).  https://arxiv.org/abs/2111.00868

The area function is described as a cascade of cylindrical tubelets,
each modeled by its electrical analog (per-length series impedance Z and
shunt admittance Y, including Fant's 1960 thermal/viscous losses and
wall-vibration admittance).  Formants are the poles of the transfer
function, found by peak-picking on the frequency response followed by a
secant (Newton) refinement in the complex frequency plane.

Conventions of the original code are kept:
  * areas/lengths are ordered glottis -> lips,
  * lengths in cm, areas in cm^2, frequencies in Hz,
  * ``no_vibration = 1`` removes the wall-vibration admittance
    ("lossless" walls setting used for the generic model and the DRM;
    the Fant model uses ``no_vibration = 0``),
  * ``TLM`` returns the two first formants, ``TLM3spec`` also returns F3
    and the 500-point magnitude spectrum used by the LPC resynthesis.
"""

import numpy as np

# Physical constants, "valeurs des constantes Mrayati" (Badin & Fant 1984)
C_AIR = 35100.0        # sound speed, cm/s
RHO_AIR = 1.14e-3      # air density, g/cm^3
LAMBDA = 5.5e-5        # thermal conductivity
ETA = 1.4              # ratio of specific heats
MU = 1.86e-4           # viscosity
CP = 0.24              # specific heat
BP = 1600.0            # wall damping, dynes.s/cm
MP = 1.4               # wall mass, g/cm^2


def spectrelec(w, A, zr, l, no_vibr_paroi=0):
    """Transfer function of the tube cascade at pulsations ``w``.

    Parameters
    ----------
    w : complex array -- pulsations 2*pi*f (may be complex for pole search)
    A : array (ntubes,) -- section areas, glottis -> lips
    zr : complex array or scalar -- radiation impedance at the lips
    l : array (ntubes,) -- section lengths
    no_vibr_paroi : bool/int -- 1 removes the wall-vibration admittance

    Returns
    -------
    H, YE : transfer function (pressure at lips / volume velocity at
    glottis, up to the radiation load) and input admittance.
    """
    # MATLAB convention: A, l as columns (ntubes,1), w as row (nbfreq,):
    # all per-tubelet quantities are outer products over the two axes.
    A = np.asarray(A).reshape(-1, 1)
    l = np.asarray(l).reshape(-1, 1)
    w = np.atleast_1d(np.asarray(w)).ravel()

    S = 2 * np.sqrt(A * np.pi)                      # perimeter of tubelet
    Lg = RHO_AIR / A * l                            # per-length inductance
    Cg = A * l / (RHO_AIR * C_AIR ** 2)             # per-length capacitance

    # Fant (1960) losses
    R_coef = np.sqrt(RHO_AIR * MU / (2 * w))
    G_coef = (ETA - 1) / (RHO_AIR * C_AIR ** 2) * np.sqrt(LAMBDA * w / (2 * CP * RHO_AIR))
    R = S * l / (A ** 2) * R_coef
    G = S * l * G_coef
    # wall vibration + viscosity along the walls
    YP_coef = 1.0 / (BP ** 2 + MP ** 2 * w ** 2)
    YP = S * l * ((BP - 1j * MP * w) * YP_coef)
    if no_vibr_paroi:
        YP = 0.0

    Z = R + 1j * Lg * w
    Y = G + 1j * Cg * w + YP

    # chain (ABCD) matrix of each tubelet, then product over the cascade
    aa = 1 + (Z * Y / 2)
    bb = -(Z + Z ** 2 * Y / 4)
    cc = -Y
    dd = aa

    aaa, bbb, ccc, ddd = aa[0], bb[0], cc[0], dd[0]
    for ind in range(len(A) - 1):
        proda = aa[ind + 1] * aaa + bb[ind + 1] * ccc
        prodb = aa[ind + 1] * bbb + bb[ind + 1] * ddd
        prodc = cc[ind + 1] * aaa + dd[ind + 1] * ccc
        prodd = cc[ind + 1] * bbb + dd[ind + 1] * ddd
        aaa, bbb, ccc, ddd = proda, prodb, prodc, prodd

    if np.isscalar(zr) and np.isinf(zr):
        H = np.zeros_like(aaa)
        YE = -ccc / ddd
    else:
        H = 1.0 / (aaa - ccc * zr)
        YE = -(aaa - ccc * zr) / (bbb - ddd * zr)
    return H, YE


def aire2spectre_oral(A, L, nbfreq, Fmax, Fmin, no_vibration):
    """Oral-tract transfer function sampled on the frequency grid.

    Returns H (linear scale) and the input admittance.  With ``Fmin`` and
    ``Fmax`` given as the *same complex scalar* and ``nbfreq = 1`` (as done
    during the pole search), the response is evaluated at that single
    complex frequency.
    """
    if Fmin is None:
        f = np.linspace(Fmax / nbfreq, Fmax, nbfreq)
    else:
        f = np.linspace(Fmin, Fmax, nbfreq)
    w = 2 * np.pi * f          # complex f allowed (complex-frequency pole search)

    # radiation impedance at the lips (last section)
    Zr = (RHO_AIR / (2 * np.pi * C_AIR) * (w ** 2)
          + 1j * 8 * RHO_AIR / (3 * np.pi * np.sqrt(np.pi * A[-1])) * w)
    return spectrelec(w, A, Zr, L, no_vibration)


def aire2spectre_cor_oral(A, L, nbfreq, Fmax, Fmin, F_form, NF, no_vibration):
    """Transfer function corrected for the ``NF`` already-found poles.

    Each pole pair (SI1, conj(SI1)) stored in the first ``NF`` rows of
    ``F_form`` (rows: [0, frequency, bandwidth]) is divided out, so that
    the peak-picking / secant search moves on to the next pole.
    """
    H_eq = aire2spectre_oral(A, L, nbfreq, Fmax, Fmin, no_vibration)[0]

    f = np.linspace(Fmin, Fmax, nbfreq)
    if np.iscomplexobj(f):
        SI = 2 * np.pi * 1j * f
    else:
        SI = 1j * 2 * np.pi * f

    F_form = np.asarray(F_form)
    for I in range(NF):
        SI1 = complex(F_form[I, 2] * np.pi, F_form[I, 1] * 2 * np.pi)
        H_eq = H_eq * ((SI - SI1) * (SI - np.conj(SI1))) / (SI1 * np.conj(SI1))
    return H_eq


def peakpick(signal):
    """Local maxima/minima of a 1-D signal (MATLAB ``peakpick``)."""
    signal = np.asarray(signal)
    left = signal[:-2]
    mid = signal[1:-1]
    right = signal[2:]
    is_max = (left < mid) & (right < mid)
    is_min = (left > mid) & (right > mid)
    ind = np.arange(1, len(signal) - 1)
    ind_maxi = ind[is_max]
    ind_mini = ind[is_min]
    if ind_maxi.size == 0:
        ind_maxi = np.array([int(np.argmax(signal))])
    if ind_mini.size == 0:
        ind_mini = np.array([int(np.argmin(signal))])
    return (signal[ind_maxi], ind_maxi, signal[ind_mini], ind_mini)


def nraph_oral(FE, BNPE, ITERMX, FMAX, A, L, F_form, NF, no_vibration):
    """Secant search of one pole of H in the complex frequency plane.

    Starting from the complex pulsation SI = -BNPE*pi + j*2*pi*FE, two
    evaluations of the pole-corrected inverse response are used to take a
    complex-secant step; convergence when the step is < SEUIL (0.3).
    Returns (F, BP): frequency and bandwidth of the pole, or (nan, nan).
    """
    DELTAS = complex(30, 30)
    SEUIL = 0.3
    SI = complex(-BNPE * np.pi, FE * 2 * np.pi)
    F, BP = np.nan, np.nan
    ITER = 1
    while ITER < ITERMX:
        with np.errstate(all='ignore'):
            FIcx = -1j * SI / (2 * np.pi)
            Q = 1.0 / complex(aire2spectre_cor_oral(A, L, 1, FIcx, FIcx, F_form, NF, no_vibration)[0])
            SIP = SI + DELTAS
            FIPcx = -1j * SIP / (2 * np.pi)
            QP = 1.0 / complex(aire2spectre_cor_oral(A, L, 1, FIPcx, FIPcx, F_form, NF, no_vibration)[0])
            Q1D = (QP - Q) / DELTAS
            SISU = SI - (Q / Q1D)
        if abs(SISU - SI) < SEUIL:
            F = np.imag(SISU) / (2 * np.pi)
            BP = -np.real(SISU) / np.pi
            return float(F), float(BP)
        SI = complex(SISU)
        ITER += 1
    return F, BP


def vtn2frm_ftr_oral(A, L, nbfreq, Fmax, Fmin, no_vibration):
    """Find the poles (formants) of the transfer function.

    Peak-picking on the (pole-corrected) response gives starting points;
    ``nraph_oral`` refines them in the complex plane.  The original
    two-pass structure of the MATLAB code is kept: first a sweep from the
    first peak upward, then any first-pass peaks not yet accounted for.
    Returns (H_eq, F_form1): plain spectrum and first three formant
    frequencies (padded with zeros), as the compiled ``vtn2frm_ftr_oral``.
    """
    f = np.linspace(Fmin, Fmax, nbfreq)
    FE = 150.0
    BNPE = 50.0
    FINC = 100.0
    ITERMX = 100
    seuil2 = 10.0
    FMAX = Fmax

    F_form = np.zeros((0, 3))
    H_eq = aire2spectre_cor_oral(A, L, nbfreq, Fmax, Fmin, F_form, 0, no_vibration)
    _, ind_maxi, _, _ = peakpick(20 * np.log10(np.abs(H_eq)))
    frq_max_first = f[ind_maxi]

    F = 0.0
    NF = 0
    while F <= FMAX and NF < 40:                    # NF guard for safety
        H_eq = aire2spectre_cor_oral(A, L, nbfreq, Fmax, Fmin, F_form, NF, no_vibration)
        _, ind_maxi, _, _ = peakpick(20 * np.log10(np.abs(H_eq)))
        FE = f[ind_maxi[0]]                         # first peak of the corrected response
        F, BP = nraph_oral(FE, BNPE, ITERMX, FMAX, A, L, F_form, NF, no_vibration)
        NF += 1
        F_form = np.vstack([F_form, [0.0, F, BP]])
        if not np.isfinite(F):
            break

    # second pass on the peaks of the initial (uncorrected) response
    seuil3 = 20.0
    keep = np.abs(frq_max_first[:, None] - F_form[None, :, 1]) >= seuil3
    frq_max = frq_max_first[keep.all(axis=1)]
    for FE in frq_max:
        F, BP = nraph_oral(FE, BNPE, ITERMX, FMAX, A, L, F_form, NF, no_vibration)
        NF += 1
        F_form = np.vstack([F_form, [0.0, F, BP]])

    # clean up (compiled version): out-of-range rows flagged at -1, sorted,
    # valid poles pulled to the top, redundant poles (< seuil2 apart) merged.
    # MATLAB keeps row i when F(i+1)-F(i) >= seuil2 (i = last of its cluster)
    # plus the final row:  [ind_F_diffOK ; ind_F_diffOK(end)+1].
    F_form = F_form[:NF]
    bad = ~np.isfinite(F_form[:, 1]) | (F_form[:, 1] < 0) | (F_form[:, 1] > FMAX)
    F_form[bad, 1] = -1.0
    F_form = F_form[np.argsort(F_form[:, 1])]
    valid = (F_form[:, 1] >= 0) & (F_form[:, 1] <= Fmax)
    F_form = F_form[valid]
    if F_form.shape[0] > 1:
        d = np.diff(F_form[:, 1])
        ok0 = np.where(d >= seuil2)[0]              # 0-based diff indices
        kept1 = np.concatenate([ok0 + 1, [ok0[-1] + 2]])   # 1-based rows
        F_form = F_form[kept1 - 1]
    else:
        F_form = F_form[:1]

    # final output: PLAIN (uncorrected) spectrum + first three formants
    H_eq = aire2spectre_oral(A, L, nbfreq, Fmax, Fmin, no_vibration)[0]
    F_form1 = np.zeros(3)
    F_form1[:min(3, F_form.shape[0])] = F_form[:3, 1]
    return H_eq, F_form1


def TLM(area, no_vibration):
    """First two formants of an area function -- MATLAB ``TLM``.

    Parameters
    ----------
    area : (ntubes, 2) array -- columns: [length (cm), area (cm^2)],
           ordered glottis -> lips.
    no_vibration : 1 for lossless walls, 0 for wall vibration.
    """
    L, A = np.asarray(area)[:, 0], np.asarray(area)[:, 1]
    nbfreq, Fmax = 500, 10000.0
    Fmin = Fmax / (nbfreq - 1)
    _, F_form1 = vtn2frm_ftr_oral(A, L, nbfreq, Fmax, Fmin, no_vibration)
    return float(F_form1[0]), float(F_form1[1])


def TLM3spec(area, no_vibration):
    """First three formants + magnitude spectrum -- MATLAB ``TLM3spec``.

    The 500-point spectrum (linear scale, 20 Hz resolution between 20 Hz
    and 10 kHz) is the input of the LPC vowel resynthesis (``voysynth``).
    """
    L, A = np.asarray(area)[:, 0], np.asarray(area)[:, 1]
    nbfreq, Fmax = 500, 10000.0
    Fmin = Fmax / (nbfreq - 1)
    H_eq, F_form1 = vtn2frm_ftr_oral(A, L, nbfreq, Fmax, Fmin, no_vibration)
    return float(F_form1[0]), float(F_form1[1]), float(F_form1[2]), np.abs(H_eq)
