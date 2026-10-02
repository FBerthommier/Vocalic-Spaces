"""Figure 4 -- (a) design of Fant's model (redrawn from Badin et al. 1990,
Fig. 1) with the geometrical parameters of Table III; (b) C1 simulation
together with the C2 vowel space at rho = 1 (coordination function).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from figures_common import load_or_compute, load_mc, savefig, IPA
from simu_fant import fant_cycle
from vsmodel import models as mdl

# Table III constants (n = 200 tubelets)
N = 200
LC = 0.9 * N
LIP = N - LC
LT = 0.3 * N


def panel_a(ax):
    """Schematic of Fant's model (redrawn from Badin et al. 1990, Fig. 1).

    Flat upper wall; the airway narrows at the constriction (lower wall
    raised, channel of area Ac below the upper wall); reduced lip opening
    Al at the right end.  Section boundaries and the Xc arrow with the
    dashed line at the constriction center, as in the submitted figure.
    """
    L1 = 60.0                        # back cavity [0, L1]
    l = LT                           # constriction [L1, L1+l], l = 60
    xc = L1 + l / 2.0                # constriction center = Omega(Xc)
    h = 0.72                         # constriction channel lower bound
    lip_h = 0.66                     # lip opening height
    xlip = LC                        # 180

    # airway polygon (thick outline, white interior)
    poly = [(0, 0), (L1, 0), (L1, h), (L1 + l, h), (L1 + l, 0),
            (xlip, 0), (xlip, lip_h), (N, lip_h), (N, 1.0), (0, 1.0)]
    ax.add_patch(plt.Polygon(poly, closed=True, facecolor='white',
                             edgecolor='k', lw=3, zorder=2))

    # dashed line at the constriction center + Xc double arrow
    ax.plot([xc, xc], [h, -0.95], 'k--', lw=1.0)
    ax.annotate('', xy=(xc, -0.62), xytext=(0, -0.62),
                arrowprops=dict(arrowstyle='<->', color='k', lw=1.2))
    ax.text(xc / 2 + 3, -1.12, r'$X_c$', fontsize=14)

    # L and Lc double arrows above
    ax.annotate('', xy=(N, 1.62), xytext=(0, 1.62),
                arrowprops=dict(arrowstyle='<->', color='k', lw=1.2))
    ax.text(N / 2, 1.70, r'$L$', fontsize=14, ha='center')
    ax.annotate('', xy=(LC, 1.28), xytext=(0, 1.28),
                arrowprops=dict(arrowstyle='<->', color='k', lw=1.2))
    ax.text(LC / 2, 1.36, r'$L_c$', fontsize=14, ha='center')

    # section labels and boundary arrows below the tube
    ylab = -0.30
    ax.plot([0, N], [ylab, ylab], 'k-', lw=0.8)
    for x in (0, L1, L1 + l, LC, N):
        ax.plot([x, x], [ylab - 0.06, ylab + 0.06], 'k-', lw=0.8)
    secs = [(0, L1, 'Back Cavity'), (L1, L1 + l, 'Constriction'),
            (L1 + l, LC, 'Front Cavity'), (LC, N, 'Lips')]
    for a, b, lab in secs:
        ax.annotate('', xy=(b, ylab - 0.16), xytext=(a, ylab - 0.16),
                    arrowprops=dict(arrowstyle='<->', color='k', lw=1.0))
        rot = 90 if lab == 'Constriction' else 0
        ax.text((a + b) / 2, ylab - 0.52 if rot else ylab - 0.48, lab,
                fontsize=10.5, ha='center', va='top', rotation=rot)

    # glottis, Ac, Al
    ax.text(-14, 0.5, 'Glottis', fontsize=11, rotation=90,
            ha='center', va='center')
    ax.text(xc, (1.0 + h) / 2, r'$A_c$', fontsize=14, ha='center',
            va='center')
    ax.annotate('', xy=(N + 6, 0), xytext=(N + 6, lip_h),
                arrowprops=dict(arrowstyle='<->', color='k', lw=1.0))
    ax.text(N + 16, lip_h / 2, r'$A_l$', fontsize=14, ha='left',
            va='center')

    ax.set_xlim(-24, N + 26)
    ax.set_ylim(-1.35, 2.45)
    ax.axis('off')
    ax.text(N / 2, 2.25, 'FANT', fontsize=20, ha='center', va='center',
            fontweight='bold')
    ax.text(-20, 2.25, '(a)', fontsize=15, ha='center', va='center')


def main(reference=False):
    d = load_or_compute('fant_cycle_py.npz', lambda: fant_cycle())
    mc, src = load_mc('fant', reference=reference)
    indb = np.array([1, 11, 16, 21, 31, 41, 46, 51]) - 1

    fig = plt.figure(figsize=(12.5, 6.4))
    axa = fig.add_axes([0.03, 0.06, 0.44, 0.86])
    panel_a(axa)

    axb = fig.add_axes([0.55, 0.14, 0.42, 0.78])
    axb.plot(mc['f2r'], mc['f1r'], '.', color='r', ms=2, alpha=0.35, label='C1')
    axb.plot(d['f2'], d['f1'], 'b-', lw=1.5, label='C2 (TLM)')
    axb.plot(d['f2'][indb], d['f1'][indb], 'o', color='b', ms=7)
    for lab, k in zip(mdl.VOWELS, indb):
        axb.annotate(IPA[lab], (d['f2'][k], d['f1'][k]), textcoords='offset points',
                     xytext=(7, -3), fontsize=13)
    axb.plot(d['f2b'], d['f1b'], 'bo', ms=8)
    axb.annotate('(f1n, f2n)', (d['f2b'], d['f1b']), xytext=(-70, 8),
                 textcoords='offset points')
    axb.set_xlim(600, 2300)
    axb.set_ylim(225, 675)
    axb.invert_yaxis()
    axb.set_xlabel('F2 (Hz)')
    axb.set_ylabel('F1 (Hz)')
    axb.legend(loc='upper right', fontsize=9)
    axb.set_title('(b) Fant model  (C1 data: %s)' % src)
    savefig(fig, 'Figure4.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB C1 draws instead of the Python ones')
    main(ap.parse_args().reference)
