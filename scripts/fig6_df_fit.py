"""Figure 6 -- estimates of the deviations from neutral for DRM C2
(faithful to the Matlab script ``figdf.m``).

The relative formant deviations follow  d~f1 = -a1/2 - (a2/2)^2  and
d~f2 = -a2/2 : the quadratic term biases the Schroeder-Ehrenfest
relation, bringing df1 close to the bisecting line.
"""

import numpy as np
import matplotlib.pyplot as plt

from figures_common import load_or_compute, load_mc, savefig
from simu_drm import drm_cycle


def main(reference=False):
    d = load_or_compute('drm_cycle_py.npz', lambda: drm_cycle())
    mc, src = load_mc('drm', reference=reference)
    f1b, f2b = float(mc['f1b']), float(mc['f2b'])

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.grid(True, alpha=0.4)
    # bisecting line
    ax.plot([-1, 1], [-1, 1], 'k--', lw=0.8)
    # C2 cloud (random rho, theta):  df1 vs estimate with the quadratic bias
    ax.plot((mc['f1s'] - f1b) / f1b, -(mc['a1s'] / 2) - (mc['a2s'] / 2) ** 2,
            '.', color='0.6', ms=2)
    ax.plot((mc['f2s'] - f2b) / f2b, -(mc['a2s'] / 2), '.', color='mistyrose',
            ms=2)
    # C2 cycle: with bias (thick) and without (thin)
    ax.plot((d['f1'] - f1b) / f1b, -(d['a1'] / 2) - (d['a2'] / 2) ** 2,
            'b-', lw=2, label=r'$df_1$ vs $-\frac{a_1}{2}-(\frac{a_2}{2})^2$')
    ax.plot((d['f1'] - f1b) / f1b, -(d['a1'] / 2), 'b-', lw=1, alpha=0.6,
            label=r'$df_1$ vs $-a_1/2$ (SE, no bias)')
    ax.plot((d['f2'] - f2b) / f2b, -(d['a2'] / 2), 'r-', lw=2,
            label=r'$df_2$ vs $-a_2/2$')
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    ax.set_aspect('equal')
    ax.set_xlabel(r'deviation from neutral  $df_i = (f_i - f_{in})/f_{in}$')
    ax.set_ylabel('estimate from $(\\tilde a_1, \\tilde a_2)$')
    ax.legend(loc='upper left', fontsize=9)
    ax.set_title('DRM C2  [data: %s]' % src)
    savefig(fig, 'Figure6.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB draws instead of the Python ones')
    main(ap.parse_args().reference)
