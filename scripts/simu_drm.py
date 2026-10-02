"""Simulation of the DRM 4-tube model driven by the coordination function
(Fig. 3 / Mm. 3-4).

Python counterpart of the supplementary Matlab script ``simuDRM.m``
(F. Berthommier, 11/24/2021): the 4 section areas follow the coordination
function P(rho, theta) = Omega + rho Psi1 cos(Psi2 - theta) with the
{Omega, Psi1, Psi2} obtained from the {i, j, k} of the generic model
(Table II); the cycle rho = 1, theta = 0 -> 2 pi gives the C2 vowel space
together with the Schroeder-Ehrenfest estimates ("SE") and the biased
estimate ("est").

Usage:  python simu_drm.py [--quick]
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from vsmodel import tlm, models as mdl   # noqa: E402

L = 17.5
N = 120             # tubelets (multiple of 6, adapted to the DRM cut)


def drm_cycle(no_vibration=1, n=N, pastheta=np.pi / 60, rho=1.0):
    """C2 cycle of the DRM + Fourier coefficients + formant estimates."""
    Om, Ps1, Ps2 = mdl.drm_setup()
    ntheta = int(round(2 * np.pi / pastheta)) + 2
    theta = np.linspace(0, pastheta * (ntheta - 1), ntheta)

    a1 = np.zeros(ntheta); a2 = np.zeros(ntheta)
    a10 = np.zeros(ntheta); a20 = np.zeros(ntheta)
    f1 = np.zeros(ntheta); f2 = np.zeros(ntheta)
    for k, th in enumerate(theta):
        Pval = Om + rho * Ps1 * np.cos(Ps2 - th)
        P = mdl.drm_P(Pval, n=n)
        a1[k], a2[k] = mdl.fourier_a12(mdl.expif(P), n=n)
        a10[k], a20[k] = mdl.fourier_a12(P, n=n)
        f1[k], f2[k] = tlm.TLM(mdl.drm_area(Pval, L=L, n=n), no_vibration)

    f1b, f2b = tlm.TLM(np.column_stack([np.full(n, L / n), np.ones(n)]), no_vibration)
    # formant estimates: "est" (quadratic bias on F1) and "SE" (linear),
    # from the rectified (a1, a2) and raw (a10, a20) coefficients
    f1est = f1b - f1b * ((a1 / 2) + (a2 / 2) ** 2)
    f2est = f2b - f2b * a2 / 2
    f1est0 = f1b - f1b * a1 / 2          # SE from (a1, a2)
    f2est0 = f2b - f2b * a2 / 2
    f1est1 = f1b - f1b * a10 / 2         # SE from (a10, a20) before rectification
    f2est1 = f2b - f2b * a20 / 2
    return dict(theta=theta, a1=a1, a2=a2, a10=a10, a20=a20, f1=f1, f2=f2,
                f1b=f1b, f2b=f2b, f1est=f1est, f2est=f2est,
                f1est0=f1est0, f2est0=f2est0, f1est1=f1est1, f2est1=f2est1,
                Omega=Om, Psi1=Ps1, Psi2=Ps2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='coarser theta step (pi/15)')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    pastheta = np.pi / 15 if args.quick else np.pi / 60

    t0 = time.time()
    res = drm_cycle(pastheta=pastheta)
    print('DRM cycle: %d configurations in %.1f s' % (len(res['theta']), time.time() - t0))
    print('neutral tube: f1b=%.2f Hz  f2b=%.2f Hz' % (res['f1b'], res['f2b']))
    print('Omega=%s' % np.round(res['Omega'], 4))
    print('Psi2 =%s rad' % np.round(res['Psi2'], 4))

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                   'data', 'drm_cycle_py.npz')
    np.savez(out, **res)
    print('saved', out)


if __name__ == '__main__':
    main()
