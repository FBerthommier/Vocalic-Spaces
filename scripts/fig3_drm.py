"""Figure 3 -- DRM: C1 simulation together with the C2 vowel space at
rho = 1 (coordination function), plus the "SE" and biased "est" formant
estimates (Schroeder-Ehrenfest relation before/after rectification).
"""

import numpy as np
import matplotlib.pyplot as plt

from figures_common import load_or_compute, load_mc, savefig, IPA
from simu_drm import drm_cycle
from vsmodel import models as mdl


def main(reference=False):
    d = load_or_compute('drm_cycle_py.npz', lambda: drm_cycle())
    mc, src = load_mc('drm', reference=reference)
    indb = mdl.INDTHETAB - 1

    fig, ax = plt.subplots(figsize=(8, 8))
    # C1: 4 uncorrelated parameters uniform in Omega +- Psi1
    ax.plot(mc['f2r'], mc['f1r'], '.', color='r', ms=2, alpha=0.35, label='C1')
    # C2 cycle with the TLM
    ax.plot(d['f2'], d['f1'], 'b-', lw=1.5, label='C2 (TLM)')
    ax.plot(d['f2'][indb], d['f1'][indb], 'o', color='b', ms=7)
    # estimates from the cosine coefficients
    ax.plot(d['f2est'], d['f1est'], 'c-', lw=1.2, label='"est" (a$_1$, a$_2$)')
    ax.plot(d['f2est'][indb], d['f1est'][indb], 'o', color='c', ms=6)
    ax.plot(d['f2est0'], d['f1est0'], 'g-', lw=1.0,
            label='SE after rectif. (a$_1$, a$_2$)')
    ax.plot(d['f2est1'], d['f1est1'], 'g--', lw=1.0,
            label='SE before rectif. (a$_{10}$, a$_{20}$)')
    ax.plot(d['f2est0'][indb], d['f1est0'][indb], 'o', color='g', ms=5)
    # neutral position
    ax.plot(d['f2b'], d['f1b'], 'bo', ms=8)
    ax.annotate('(f1n, f2n)', (d['f2b'], d['f1b']), xytext=(-75, 8),
                textcoords='offset points')
    # vowel labels on the C2 loop
    for lab, k in zip(mdl.VOWELS, indb):
        ax.annotate(IPA[lab], (d['f2'][k], d['f1'][k]), textcoords='offset points',
                    xytext=(7, -3), fontsize=13)

    ax.set_xlim(250, 2650)
    ax.set_ylim(0, 1000)
    ax.invert_yaxis()
    ax.set_xlabel('F2 (Hz)')
    ax.set_ylabel('F1 (Hz)')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title('DRM  (C1 data: %s)' % src)
    savefig(fig, 'Figure3.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB C1 draws instead of the Python ones')
    main(ap.parse_args().reference)
