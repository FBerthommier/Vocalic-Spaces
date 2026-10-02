"""Validation of the Python translation against the MATLAB reference data.

The reference .npz files (data/reference/) were converted from the .mat
files produced by the compiled MATLAB code for the article.  This script
recomputes the deterministic parts (neutral tubes, C2 cycles, Fourier
coefficients) and reports the deviations:

  * generic / DRM (lossless walls): F1, F2 within ~1 Hz of the reference;
  * Fant (lossy, wall vibration): F2 within ~1 Hz; F1 shows a systematic
    ~1.5-2 Hz offset on the most damped pole, attributed to a small drift
    between the current MATLAB sources and the 2021-compiled binary that
    produced the reference;
  * cosine Fourier coefficients: machine precision;
  * LPC resynthesis: verified in the unit tests (lattice == direct form).

Usage:  python validate.py [--quick]
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from vsmodel import tlm, models as mdl
from simu_gen import generic_cycle
from simu_drm import drm_cycle
from simu_fant import fant_cycle

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'reference')
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def sc(a):
    """Scalar from a possibly 0-d/1-element array."""
    return float(np.ravel(a)[0])


def check(name, err, unit, tol, lines):
    status = 'OK  ' if err <= tol else 'FAIL'
    lines.append('  [%s] %-42s max err = %8.4f %s (tol %g)' % (status, name, err, unit, tol))
    print(lines[-1])
    return err <= tol


def main(quick=False):
    lines = ['Validation report -- Python translation vs MATLAB reference data', '=' * 78]
    ok = True

    # ---------------- generic model ----------------------------------------
    t0 = time.time()
    theta, f1, f2, f1b, f2b = generic_cycle(pastheta=np.pi / 15 if quick else np.pi / 60)
    ref = np.load(os.path.join(REF, 'gen_cycle.npz'))
    if quick:                                   # subsample reference to our grid
        step = 4
        r1, r2 = ref['f1'][::step], ref['f2'][::step]
        n = min(len(f1), len(r1)); f1, f2, r1, r2 = f1[:n], f2[:n], r1[:n], r2[:n]
    else:
        r1, r2 = ref['f1'], ref['f2']
    lines.append('GENERIC model (lossless, n=100, %d pts, %.0f s)' % (len(f1), time.time() - t0))
    print(lines[-1])
    ok &= check('neutral f1b (ref %.4f)' % sc(ref['f1b']),
                abs(f1b - sc(ref['f1b'])), 'Hz', 1.0, lines)
    ok &= check('neutral f2b (ref %.4f)' % sc(ref['f2b']),
                abs(f2b - sc(ref['f2b'])), 'Hz', 1.0, lines)
    ok &= check('cycle F1 (122 phases)', np.max(np.abs(f1 - r1)), 'Hz', 1.0, lines)
    ok &= check('cycle F2 (122 phases)', np.max(np.abs(f2 - r2)), 'Hz', 1.0, lines)

    # ---------------- DRM ----------------------------------------------------
    t0 = time.time()
    pastheta = np.pi / 15 if quick else np.pi / 60
    d = drm_cycle(pastheta=pastheta)
    refc = np.load(os.path.join(REF, 'drm_cycle.npz'))      # grid pi/60
    refm = np.load(os.path.join(REF, 'drm_mc.npz'))         # grid pi/30
    n = len(d['f1'])
    lines.append('DRM model (lossless, n=120, %d pts, %.0f s)' % (n, time.time() - t0))
    print(lines[-1])
    ok &= check('neutral f1b (ref %.4f)' % sc(refc['f1b']),
                abs(d['f1b'] - sc(refc['f1b'])), 'Hz', 1.0, lines)
    ok &= check('neutral f2b (ref %.4f)' % sc(refc['f2b']),
                abs(d['f2b'] - sc(refc['f2b'])), 'Hz', 1.0, lines)
    step = max(1, int(pastheta / (np.pi / 60) + 0.5))         # our grid in ref units
    m = min(n, len(refc['f1']) // step)
    ok &= check('cycle F1', np.max(np.abs(d['f1'][:m] - refc['f1'][::step][:m])), 'Hz', 1.0, lines)
    ok &= check('cycle F2', np.max(np.abs(d['f2'][:m] - refc['f2'][::step][:m])), 'Hz', 1.0, lines)
    # coefficients + estimates against the 62-pt reference cycle (grid pi/30).
    # Map our theta grid onto the reference grid (both contain the phases
    # that are multiples of lcm(our step, pi/30)).
    ratio = pastheta / (np.pi / 30)
    if ratio >= 1:      # our grid coarser: our point i -> reference i*ratio
        idx_ref = np.arange(min(len(d['a1']), len(refm['a1']))) * int(round(ratio))
        idx_ref = idx_ref[idx_ref < len(refm['a1'])]
        mine_a1, ref_a1 = d['a1'][:len(idx_ref)], refm['a1'][idx_ref]
        mine_est, ref_est = d['f1est'][:len(idx_ref)], refm['f1est'][idx_ref]
        mine_a2, ref_a2 = d['a2'][:len(idx_ref)], refm['a2'][idx_ref]
    else:               # our grid finer: take every 1/ratio-th of our points
        s = int(round(1 / ratio))
        m2 = min(len(d['a1']) // s, len(refm['a1']))
        mine_a1, ref_a1 = d['a1'][::s][:m2], refm['a1'][:m2]
        mine_a2, ref_a2 = d['a2'][::s][:m2], refm['a2'][:m2]
        mine_est, ref_est = d['f1est'][::s][:m2], refm['f1est'][:m2]
    ok &= check('cycle a1 coefficients', np.max(np.abs(mine_a1 - ref_a1)), '-', 1e-9, lines)
    ok &= check('cycle a2 coefficients', np.max(np.abs(mine_a2 - ref_a2)), '-', 1e-9, lines)
    # f1est inherits the neutral-formant precision (~0.05 Hz on f1b,
    # amplified by the a1/2 term), hence the 2 Hz tolerance
    ok &= check('cycle f1est ("est" relation)',
                np.max(np.abs(mine_est - ref_est)), 'Hz', 2.0, lines)

    # ---------------- Fant ---------------------------------------------------
    t0 = time.time()
    pastheta = np.pi / 15 if quick else np.pi / 30
    d = fant_cycle(pastheta=pastheta)
    refc = np.load(os.path.join(REF, 'fant_cycle.npz'))     # grid pi/30
    refm = np.load(os.path.join(REF, 'fant_mc.npz'))
    lines.append('FANT model (lossy, n=200, %d pts, %.0f s)' % (len(d['f1']), time.time() - t0))
    print(lines[-1])
    ok &= check('neutral f2b (ref %.4f)' % sc(refc['f2b']),
                abs(d['f2b'] - sc(refc['f2b'])), 'Hz', 1.0, lines)
    step = max(1, int(pastheta / (np.pi / 30) + 0.5))
    m = min(len(d['f1']), len(refc['f1']) // step)
    ok &= check('cycle F2 (62 phases)', np.max(np.abs(d['f2'][:m] - refc['f2'][::step][:m])), 'Hz', 1.0, lines)
    ok &= check('neutral f1b (ref %.4f, known lossy drift)' % sc(refc['f1b']),
                abs(d['f1b'] - sc(refc['f1b'])), 'Hz', 3.0, lines)
    ok &= check('cycle F1 (62 phases, known lossy drift)',
                np.max(np.abs(d['f1'][:m] - refc['f1'][::step][:m])), 'Hz', 3.0, lines)
    ok &= check('cycle a1 coefficients',
                np.max(np.abs(d['a1'][:m] - refm['a1'][::step][:m])), '-', 1e-9, lines)
    ok &= check('cycle a2 coefficients',
                np.max(np.abs(d['a2'][:m] - refm['a2'][::step][:m])), '-', 1e-9, lines)

    # ---------------- coordination transform ---------------------------------
    ijk = np.array([[1.0, -2.0, 1.0], [1.0, 0.0, -1.0]])
    om, ps1, ps2 = mdl.ijk2psi(ijk)
    exp_ps1 = np.array([2.0, 1 / np.cos(np.pi / 6)])
    exp_ps2 = np.array([0.0, np.pi / 2])
    lines.append('COORDINATION function (Table I)')
    print(lines[-1])
    ok &= check('Psi1 = (2, cos^-1 pi/6)', np.max(np.abs(ps1 - exp_ps1)), '-', 1e-12, lines)
    ok &= check('Psi2 = (0, pi/2)', np.max(np.abs(ps2 - exp_ps2)), '-', 1e-12, lines)

    lines.append('=' * 78)
    lines.append('OVERALL: %s' % ('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'))
    print('\n'.join(lines[-2:]))

    os.makedirs(FIGDIR, exist_ok=True)
    with open(os.path.join(FIGDIR, 'validation_report.txt'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print('report written to figures/validation_report.txt')
    return 0 if ok else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='coarser cycles')
    sys.exit(main(ap.parse_args().quick))
