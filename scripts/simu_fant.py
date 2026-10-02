"""Simulation of FANT's 4-parameter model driven by the coordination
function (Fig. 4b / Mm. 5-6).

Python counterpart of the supplementary Matlab script ``simuFANT.m``
(F. Berthommier, 11/24/2021): the four parameters {Xc, Ac, Al, L} follow
P(rho, theta) = Omega + rho Psi1 cos(Psi2 - theta) (Table III); the
n-tube reconstruction is done by rounding, merging, then soft
rectification, and the transmission line model is used WITH wall
vibration losses (lossy).

Usage:  python simu_fant.py [--quick]
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from vsmodel import tlm, models as mdl   # noqa: E402

N = 200             # tubelets (high number, avoids rounding effects)


def fant_cycle(no_vibration=0, n=N, pastheta=np.pi / 30, rho=1.0):
    """C2 cycle of Fant's model (Table III parameters)."""
    Om, Ps1, Ps2, Lc, l, Lip, A = mdl.fant_setup(n)
    ntheta = int(round(2 * np.pi / pastheta)) + 2
    theta = np.linspace(0, pastheta * (ntheta - 1), ntheta)

    a1 = np.zeros(ntheta); a2 = np.zeros(ntheta)
    a10 = np.zeros(ntheta); a20 = np.zeros(ntheta)
    f1 = np.zeros(ntheta); f2 = np.zeros(ntheta)
    Pval = Om.copy()
    L1 = int(mdl.mround(Om[0] - l / 2))
    L2 = int(Lc - l - L1)
    P = np.concatenate([A * np.ones(L1), Om[1] * np.ones(int(l)),
                        A * np.ones(L2), Om[2] * np.ones(int(Lip))])
    f1b, f2b = tlm.TLM(np.column_stack([np.full(n, Om[3] / n), mdl.expif(P)]), no_vibration)

    for k, th in enumerate(theta):
        Pval = Om + rho * Ps1 * np.cos(Ps2 - th)
        Praw = mdl.fant_P_raw(Pval, Lc, l, Lip, A)
        area = mdl.fant_area(Pval, Lc, l, Lip, A, n=n)
        f1[k], f2[k] = tlm.TLM(area, no_vibration)
        a1[k], a2[k] = mdl.fourier_a12(area[:, 1], n=n)
        a10[k], a20[k] = mdl.fourier_a12(Praw, n=n)
    return dict(theta=theta, a1=a1, a2=a2, a10=a10, a20=a20, f1=f1, f2=f2,
                f1b=f1b, f2b=f2b, Omega=Om, Psi1=Ps1, Psi2=Ps2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='coarser theta step (pi/12)')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    pastheta = np.pi / 12 if args.quick else np.pi / 30

    t0 = time.time()
    res = fant_cycle(pastheta=pastheta)
    print('FANT cycle: %d configurations in %.1f s' % (len(res['theta']), time.time() - t0))
    print('neutral config: f1b=%.2f Hz  f2b=%.2f Hz' % (res['f1b'], res['f2b']))

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                   'data', 'fant_cycle_py.npz')
    np.savez(out, **res)
    print('saved', out)


if __name__ == '__main__':
    main()
