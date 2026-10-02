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
    """Schematic of Fant's model: sliding constriction of length l."""
    L1 = 75.0                                   # Omega(Xc) - l/2 at neutral
    l = LT
    L2 = LC - l - L1
    A, Ac, Al = 1.0, -1.5, 0.5                  # raw values (before expif)
    xs = [0, L1, L1, L1 + l, L1 + l, N, N]
    ys = [A, A, Ac, Ac, A, A, Al]
    # mirrored area profile (upper/lower walls)
    ax.fill_between(xs, [y / 2 for y in ys], [-y / 2 for y in ys],
                    color='#e8e8e8', zorder=1)
    ax.plot(xs, [y / 2 for y in ys], 'k-', lw=1.5)
    ax.plot(xs, [-y / 2 for y in ys], 'k-', lw=1.5)
    ax.plot([0, N], [0, 0], 'k:', lw=0.8)
    # annotations
    arrow = dict(arrowstyle='<->', color='b', lw=1.0)
    ax.add_patch(FancyArrowPatch((0, 1.35), (L1 + l / 2, 1.35), **arrow))
    ax.text(L1 / 3, 1.5, r'$X_c$', color='b', fontsize=13)
    ax.add_patch(FancyArrowPatch((L1, -1.05), (L1 + l, -1.05), **arrow))
    ax.text(L1 + l / 2 - 8, -1.45,(r'$l$'), color='b', fontsize=13)
    ax.add_patch(FancyArrowPatch((0, -1.55), (LC, -1.55), **arrow))
    ax.text(LC / 2 - 12, -1.95, r'$L_c = 0.9\,L$', color='b', fontsize=13)
    ax.add_patch(FancyArrowPatch((0, 1.85), (N, 1.85), **arrow))
    ax.text(N / 2 - 10, 2.0, r'$L$ (through $P_4$)', color='b', fontsize=13)
    ax.annotate(r'$A_c$ (through $P_2$)', (L1 + l / 2, Ac / 2),
                xytext=(L1 + l / 2 + 30, -0.75), fontsize=12,
                arrowprops=dict(arrowstyle='->', color='b'))
    ax.annotate(r'$A_l$ (through $P_3$)', (N - 6, Al / 2),
                xytext=(N - 75, 0.95), fontsize=12,
                arrowprops=dict(arrowstyle='->', color='b'))
    ax.text(4, 0.72, 'glottis', fontsize=10)
    ax.text(N - 55, 0.72, 'lips', fontsize=10)
    ax.set_xlim(-8, N + 8)
    ax.set_ylim(-2.2, 2.3)
    ax.axis('off')
    ax.set_title('(a) Fant model design (Badin et al., 1990)', fontsize=11)


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
