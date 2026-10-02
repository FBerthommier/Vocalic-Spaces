"""Simulation of the GENERIC vocal-tract model (Fig. 1 / Mm. 1-2).

Python counterpart of the supplementary Matlab script ``simuGEN.m``
(F. Berthommier, 11/24/2021): sweeps the three-phase mixing function
theta = 0 -> 2 pi at rho = 1, computes (F1, F2) with the transmission-line
model for each configuration and saves the vowel-space cycle.

Usage:  python simu_gen.py [--quick]
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from vsmodel import tlm, models as mdl   # noqa: E402

L = 17.5            # tract length (cm)
N = 100             # tubelets


def generic_cycle(no_vibration=1, n=N, pastheta=np.pi / 60):
    """Vowel-space cycle of the generic model at rho = 1 (Eq. 3b + Eq. 4ab)."""
    ntheta = int(round(2 * np.pi / pastheta)) + 2     # 0 : pastheta : 2*pi+pastheta
    theta = np.linspace(0, pastheta * (ntheta - 1), ntheta)
    f1 = np.zeros(ntheta)
    f2 = np.zeros(ntheta)
    for k, th in enumerate(theta):
        P = mdl.generic_P(th, rho=1.0, n=n)
        f1[k], f2[k] = tlm.TLM(mdl.generic_area(P, L=L, n=n), no_vibration)
    f1b, f2b = tlm.TLM(np.column_stack([np.full(n, L / n), np.ones(n)]), no_vibration)
    return theta, f1, f2, f1b, f2b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='coarser theta step (pi/15)')
    ap.add_argument('--out', default=None, help='output .npz (default data/gen_cycle_py.npz)')
    args = ap.parse_args()
    pastheta = np.pi / 15 if args.quick else np.pi / 60

    t0 = time.time()
    theta, f1, f2, f1b, f2b = generic_cycle(pastheta=pastheta)
    print('generic cycle: %d configurations in %.1f s' % (len(theta), time.time() - t0))
    print('neutral tube: f1b=%.2f Hz  f2b=%.2f Hz' % (f1b, f2b))
    indb = mdl.INDTHETAB - 1
    for lab, k in zip(mdl.VOWELS, indb):
        print('  [%s]  theta=%5.2f pi   F1=%6.1f  F2=%7.1f'
              % (lab, theta[k] / np.pi, f1[k], f2[k]))

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                   'data', 'gen_cycle_py.npz')
    np.savez(out, theta=theta, f1=f1, f2=f2, f1b=f1b, f2b=f2b)
    print('saved', out)


if __name__ == '__main__':
    main()
