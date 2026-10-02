"""Monte-Carlo simulations of the two 4-tube models in conditions C1 and C2
(Figs. 3, 4b, 5, 6, 7).

Python counterpart of the Matlab scripts ``simufonc5randc.m`` (DRM) and
``simufonc6rand*.m`` (Fant):

* C1: the 4 parameters are uncorrelated, uniformly drawn in Omega +- Psi1;
* C2: (rho, theta) uniformly drawn, parameters set by the coordination
  function P = Omega + rho Psi1 cos(Psi2 - theta).

For each draw the TLM formants and the cosine Fourier coefficients of the
rectified area function (a1, a2 -- Eq. 8) and of the raw one (a10, a20)
are stored.

Usage:
    python monte_carlo.py drm  [--n1 1000] [--n2 5000] [--seed 0]
    python monte_carlo.py fant [--n1 1000] [--n2 5000] [--seed 0]
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from vsmodel import tlm, models as mdl   # noqa: E402


def mc_drm(n1=1000, n2=5000, seed=0, no_vibration=1, n=120, L=17.5, verbose=True):
    Om, Ps1, Ps2 = mdl.drm_setup()
    rng = np.random.RandomState(seed)
    out = {'f1b': None, 'f2b': None}
    for cond, count in (('r', n1), ('s', n2)):
        f1 = np.zeros(count); f2 = np.zeros(count)
        a1 = np.zeros(count); a2 = np.zeros(count)
        a10 = np.zeros(count); a20 = np.zeros(count)
        Pval = np.zeros((count, 4))
        rho = np.zeros(count); theta = np.zeros(count)
        t0 = time.time()
        for k in range(count):
            if cond == 'r':                                   # C1: uncorrelated
                pv = Om + Ps1 * (2 * rng.rand(4) - 1)
            else:                                             # C2: coordination
                rho[k] = rng.rand()
                theta[k] = 2 * np.pi * rng.rand()
                pv = Om + rho[k] * Ps1 * np.cos(Ps2 - theta[k])
            Pval[k] = pv
            P = mdl.drm_P(pv, n=n)
            a1[k], a2[k] = mdl.fourier_a12(mdl.expif(P), n=n)
            a10[k], a20[k] = mdl.fourier_a12(P, n=n)
            f1[k], f2[k] = tlm.TLM(mdl.drm_area(pv, L=L, n=n), no_vibration)
            if verbose and (k + 1) % 500 == 0:
                print('  DRM C%s %d/%d (%.0f s)' % ('1' if cond == 'r' else '2',
                                                     k + 1, count, time.time() - t0))
        out['f1' + cond] = f1; out['f2' + cond] = f2
        out['a1' + cond] = a1; out['a2' + cond] = a2
        out['a10' + cond] = a10; out['a20' + cond] = a20
        out['Pval' + cond] = Pval
        if cond == 's':
            out['rhos'] = rho; out['thetas'] = theta
    out['f1b'], out['f2b'] = tlm.TLM(
        np.column_stack([np.full(n, L / n), np.ones(n)]), no_vibration)
    return out


def mc_fant(n1=1000, n2=5000, seed=0, no_vibration=0, n=200, verbose=True):
    Om, Ps1, Ps2, Lc, l, Lip, A = mdl.fant_setup(n)
    rng = np.random.RandomState(seed)
    out = {'f1b': None, 'f2b': None}
    for cond, count in (('r', n1), ('s', n2)):
        f1 = np.zeros(count); f2 = np.zeros(count)
        a1 = np.zeros(count); a2 = np.zeros(count)
        a10 = np.zeros(count); a20 = np.zeros(count)
        Pval = np.zeros((count, 4))
        rho = np.zeros(count); theta = np.zeros(count)
        t0 = time.time()
        for k in range(count):
            if cond == 'r':                                   # C1: uncorrelated
                pv = Om + Ps1 * (2 * rng.rand(4) - 1)
            else:                                             # C2: coordination
                rho[k] = rng.rand()
                theta[k] = 2 * np.pi * rng.rand()
                pv = Om + rho[k] * Ps1 * np.cos(Ps2 - theta[k])
            Pval[k] = pv
            Praw = mdl.fant_P_raw(pv, Lc, l, Lip, A)
            area = mdl.fant_area(pv, Lc, l, Lip, A, n=n)
            a1[k], a2[k] = mdl.fourier_a12(area[:, 1], n=n)
            a10[k], a20[k] = mdl.fourier_a12(Praw, n=n)
            f1[k], f2[k] = tlm.TLM(area, no_vibration)
            if verbose and (k + 1) % 500 == 0:
                print('  FANT C%s %d/%d (%.0f s)' % ('1' if cond == 'r' else '2',
                                                      k + 1, count, time.time() - t0))
        out['f1' + cond] = f1; out['f2' + cond] = f2
        out['a1' + cond] = a1; out['a2' + cond] = a2
        out['a10' + cond] = a10; out['a20' + cond] = a20
        out['Pval' + cond] = Pval
        if cond == 's':
            out['rhos'] = rho; out['thetas'] = theta
    Praw = mdl.fant_P_raw(Om, Lc, l, Lip, A)
    out['f1b'], out['f2b'] = tlm.TLM(
        np.column_stack([np.full(n, Om[3] / n), mdl.expif(Praw)]), no_vibration)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model', choices=['drm', 'fant'])
    ap.add_argument('--n1', type=int, default=1000, help='C1 draws (uncorrelated)')
    ap.add_argument('--n2', type=int, default=5000, help='C2 draws (coordination)')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    if args.model == 'drm':
        out = mc_drm(args.n1, args.n2, args.seed)
        default = 'mc_drm_py.npz'
    else:
        out = mc_fant(args.n1, args.n2, args.seed)
        default = 'mc_fant_py.npz'
    out_path = args.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'data', default)
    np.savez(out_path, **out)
    print('saved', out_path)


if __name__ == '__main__':
    main()
