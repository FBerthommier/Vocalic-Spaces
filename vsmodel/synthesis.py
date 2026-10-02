"""LPC vowel resynthesis from the TLM spectrum (``voysynth`` chain).

Translation of the ICP-Grenoble MATLAB code by Laurent Girin, adapted by
F. Berthommier (``voysynthx.m``, ``Hfreq2lpc.m``, ``gen_src_3.m``,
``f_lpc_exc2sig.m``, 2009-2021).  Used to produce the multimedia sound
files Mm. 1, 4 and 6 of the article: a 0.5 s vowel is synthesized by
filtering a glottal excitation through an order-30 LPC filter fitted on
the 500-point TLM transfer-function magnitude.

Note on fidelity: in the 2021 supplement chain, ``f_lpc_exc2sig`` loads
the fixed excitation ``excit.mat`` (10 000 samples at 20 kHz), which
overrides its ``e`` argument.  The same behaviour is kept here: pass the
path of ``excit.mat`` (shipped in ``data/``), or your own array.
"""

import numpy as np
from scipy.signal import lfilter

FS = 20000          # sampling rate (Hz)
LPC_ORDER = 30      # LPC order (a bit higher than usual, per original code)


def Hfreq2lpc(H_eq, p):
    """LPC filter A(z) fitted on a (positive-frequency) magnitude response.

    The autocorrelation sequence is the IFFT of the power spectrum
    completed by Hermitian symmetry, then Levinson-Durbin gives the
    prediction coefficients and the reflection coefficients.

    Returns (a, k, g): prediction polynomial (a[0] = 1), reflection
    coefficients, and gain.
    """
    H = np.asarray(H_eq, dtype=float)
    N = len(H)
    Hr = np.concatenate([H, [H[-1]], H[1:N - 1][::-1]])
    R = np.real(np.fft.ifft(np.abs(Hr) ** 2))[:p + 1]

    # Levinson-Durbin with reflection coefficients
    a = np.zeros(p + 1)
    a[0] = 1.0
    k = np.zeros(p)
    E = R[0]
    for i in range(1, p + 1):
        acc = R[i] + (np.dot(a[1:i], R[i - 1:0:-1]) if i > 1 else 0.0)
        lam = -acc / E
        a_prev = a.copy()
        for j in range(1, i):
            a[j] = a_prev[j] + lam * a_prev[i - j]
        a[i] = lam
        k[i - 1] = lam
        E *= (1 - lam ** 2)
    g = np.sqrt(np.dot(a, R[:p + 1]))
    return a, k, g


def gen_src_3(L, fs, F0, flag_stretch=0):
    """Glottal source with Rosenberg (1971) pulses following the F0 track."""
    F0 = np.asarray(F0, dtype=float)
    M = len(F0)
    if M != L:
        F0 = np.interp(np.arange(L) / (L - 1), np.arange(M) / (M - 1), F0)

    def fix(x):
        return int(x)                      # MATLAB fix: truncate toward zero

    if flag_stretch == 0:
        F0_max = F0.max()
        T = fix(fs / F0_max / 0.8)
        Tp = fix(T / 3)
        Tn = fix(T / 4)
        t1 = np.arange(1, Tp + 1)
        t2 = np.arange(Tp + 1, Tp + Tn + 1)
        pulse = np.concatenate([
            3 * (t1 / Tp) ** 2 - 2 * (t1 / Tp) ** 3,
            1 - ((t2 - Tp) / Tn) ** 2])
        train = np.zeros(L)
        ind = 0
        while True:
            train[ind] = 1
            ind = ind + fix(fs / F0[ind])
            if ind >= L:
                break
        s = np.convolve(train, pulse)[:L]
    else:
        s = np.zeros(L)
        ind = 0
        while True:
            T = fix(fs / F0[ind])
            Tp = fix(T / 3)
            Tn = fix(T / 6)
            t1 = np.arange(1, Tp + 1)
            t2 = np.arange(Tp + 1, Tp + Tn + 1)
            pulse = np.concatenate([
                3 * (t1 / Tp) ** 2 - 2 * (t1 / Tp) ** 3,
                1 - ((t2 - Tp) / Tn) ** 2])
            s[ind:ind + Tp + Tn] = pulse
            ind = ind + fix(fs / F0[ind])
            if ind >= L:
                break
    return s


def f_lpc_exc2sig(e, step_size, trans_size, fs, LPC, excit=None):
    """Signal synthesis through the time-varying lattice LPC filter.

    ``LPC`` rows are prediction polynomials (one per frame of
    step_size seconds); the reflection coefficients inter-frame transition
    lasts trans_size seconds.  If ``excit`` is given (array), it replaces
    ``e`` -- this mirrors the original code loading ``excit.mat``.
    """
    if excit is not None:
        e = np.asarray(excit, dtype=float)
    N_step = int(step_size * fs)
    N_trans = int(trans_size * fs)
    L = len(e)
    M = L // N_step
    LPC = np.atleast_2d(LPC)
    if M != LPC.shape[0]:
        raise ValueError('LPC frame count %d != %d' % (LPC.shape[0], M))

    p = LPC.shape[1] - 1
    sig = np.zeros(L)
    e_f = np.zeros(p + 2)                  # forward lattice buffer
    e_b = np.zeros(p + 2)                  # backward lattice buffer

    def poly2rc(a):
        """Reflection coefficients from prediction polynomial (poly2rc).

        Step-down:  a_{j-1}[i] = (a_j[i] - k_j a_j[j-i]) / (1 - k_j^2),
        k_j = a_j[j],  i = 0..j-1.
        """
        a = np.array(a, dtype=float)
        pp = len(a) - 1
        kk = np.zeros(pp)
        for j in range(pp, 0, -1):
            kk[j - 1] = a[j]
            a = (a[:j] - kk[j - 1] * a[1:j + 1][::-1]) / (1 - kk[j - 1] ** 2)
        return kk

    k_P = poly2rc(LPC[0])
    for m in range(M):
        k_C = poly2rc(LPC[m])
        dk = (k_C - k_P) / N_trans
        for ind in range(1, N_step + 1):
            if ind < N_trans + 1:
                k_P = k_P + dk
            # NOTE: Python is 0-based; buffers indexed as in MATLAB with +1.
            e_f[p + 1] = e[ind - 1 + m * N_step]
            for j in range(p, 0, -1):
                e_f[j] = e_f[j + 1] - k_P[j - 1] * e_b[j]
                e_b[j + 1] = e_b[j] + k_P[j - 1] * e_f[j]
            sig[ind - 1 + m * N_step] = e_f[1]
            e_b[1] = e_f[1]
    return sig


def voysynth(spec, excit=None, fs=FS):
    """Vowel signal from a TLM spectrum (MATLAB ``voysynthx``).

    ``spec``: 500-point TLM magnitude spectrum.  ``excit``: optional
    replacement excitation (the original chain loads ``excit.mat``,
    overriding the generated source).  Returns the normalized signal.
    """
    M = 100                                  # number of spectra
    p = LPC_ORDER
    T_s = 5e-3                               # frame period (s)
    L_t = M * T_s                            # 0.5 s
    N_t = int(L_t * fs)

    a, k, g = Hfreq2lpc(spec, p)

    # F0 trajectory: spline through (0, 100), (N/3, 125), (N, 80) + jitter
    from scipy.interpolate import CubicSpline
    xs = np.array([0, N_t // 3, N_t - 1])
    F0 = CubicSpline(xs, [100.0, 125.0, 80.0])(np.arange(N_t))
    rng = np.random.RandomState(0)
    F0 = F0 + rng.randn(N_t) * F0.mean() / 100 * 0.5

    e = gen_src_3(N_t, fs, F0, 0)
    LPC = np.tile(a, (M, 1))
    sig = f_lpc_exc2sig(e, T_s, 1e-3, fs, LPC, excit=excit)

    han = np.hanning(len(sig))               # MATLAB hanning (symmetric)
    sig = lfilter([1, -0.9375], [1], han * sig)
    return sig / (1.01 * np.max(np.abs(sig)))


def wavwrite(sig, fs, path):
    """16-bit PCM wav (like MATLAB wavwrite on double data)."""
    from scipy.io import wavfile
    sig = np.asarray(sig, dtype=float)
    sig = np.clip(sig, -1, 1)
    wavfile.write(path, fs, np.int16(np.round(sig * 32767)))
