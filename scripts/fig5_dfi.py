"""Figure 5 -- deviations from the neutral position as a function of the
cosine Fourier coefficients, in conditions C1 and C2: (a) DRM, (b) Fant.

Abscissa: Fourier coefficients (a1, a2) of the rectified area function.
Ordinate: relative deviations dfi = (fi - fin) / fin.  C2 configurations
(gray circles) collapse on surfaces whereas C1 points (dots) scatter:
the many-to-one relation disappears under the coordination function.
"""

import numpy as np
import matplotlib.pyplot as plt

from figures_common import load_mc, savefig


def panel(ax, mc, f1b, f2b, title):
    df1r = (mc['f1r'] - f1b) / f1b
    df2r = (mc['f2r'] - f2b) / f2b
    df1s = (mc['f1s'] - f1b) / f1b
    df2s = (mc['f2s'] - f2b) / f2b
    # C2 "surfaces"
    ax.plot(mc['a1s'], df1s, 'o', color='0.62', ms=2.0, alpha=0.5,
            label=r'C2:  $\tilde a_1$ vs $df_1$')
    ax.plot(mc['a2s'], df2s, 'o', color='0.62', ms=2.0, alpha=0.5,
            label=r'C2:  $\tilde a_2$ vs $df_2$')
    # C1 dots
    ax.plot(mc['a1r'], df1r, '.', color='b', ms=3, alpha=0.6,
            label=r'C1:  $\tilde a_1$ vs $df_1$')
    ax.plot(mc['a2r'], df2r, '.', color='r', ms=3, alpha=0.6,
            label=r'C1:  $\tilde a_2$ vs $df_2$')
    ax.axhline(0, color='0.7', lw=0.6)
    ax.axvline(0, color='0.7', lw=0.6)
    ax.set_xlabel(r'cosine Fourier coefficients  $\tilde a_1$,  $\tilde a_2$')
    ax.set_ylabel(r'$df_1$,  $df_2$')
    ax.legend(loc='lower left', fontsize=8)
    ax.set_title(title)


def main(reference=False):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    src = ''
    for ax, model, title in ((axes[0], 'drm', '(a) DRM'),
                             (axes[1], 'fant', '(b) Fant model')):
        mc, src = load_mc(model, reference=reference)
        panel(ax, mc, float(mc['f1b']), float(mc['f2b']), title)
    fig.suptitle('dfi in conditions C1 (dots) and C2 (gray surfaces)'
                 '   [data: %s]' % src, y=1.0)
    savefig(fig, 'Figure5.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB draws instead of the Python ones')
    main(ap.parse_args().reference)
