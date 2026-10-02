"""Figure 1 -- vowel space of the GENERIC model and its 8 characteristic
vowels, with the area function of each (pairs (a1, a2) of Eq. 3b).

Center: the (F2, F1) cycle at rho = 1 with the 8 vowels and the neutral
position.  Ring panels: the rectified area functions A(x) of the 8 vowels
placed at their angular position in the vowel space.
"""

import numpy as np
import matplotlib.pyplot as plt

from figures_common import figdir, load_or_compute, savefig, IPA
from simu_gen import generic_cycle
from vsmodel import models as mdl

L = 17.5
N = 100
A_PAIRS = ['(2, 0)', '(1, 1)', '(0, α)', '(−1, 1)', '(−2, 0)',
           '(−1, −1)', '(0, −α)', '(1, −1)']


def main():
    d = load_or_compute('gen_cycle_py.npz',
                        lambda: dict(zip(('theta', 'f1', 'f2', 'f1b', 'f2b'),
                                         generic_cycle())))
    f1, f2 = d['f1'], d['f2']
    indb = mdl.INDTHETAB - 1
    thetas_char = 2 * np.pi * np.array([0, 1 / 3, 1 / 2, 2 / 3, 1, 4 / 3, 3 / 2, 5 / 3]) / 2

    fig = plt.figure(figsize=(12.5, 10.5))

    # --- center: vowel space -------------------------------------------------
    axc = fig.add_axes([0.36, 0.36, 0.30, 0.30])
    axc.plot(f2, f1, 'b-', lw=1.5)
    axc.plot(f2[indb], f1[indb], 'o', color='b', ms=7)
    for lab, k in zip(mdl.VOWELS, indb):
        axc.annotate(IPA[lab], (f2[k], f1[k]), textcoords='offset points',
                     xytext=(6, -2), fontsize=13)
    axc.plot(d['f2b'], d['f1b'], 'bo', ms=7)
    axc.annotate('(f1n, f2n)', (d['f2b'], d['f1b']), textcoords='offset points',
                 xytext=(-70, 8), fontsize=9)
    axc.set_xlim(600, 2400)
    axc.set_ylim(150, 1000)
    axc.invert_yaxis()
    axc.set_xlabel('F2 (Hz)', fontsize=9)
    axc.set_ylabel('F1 (Hz)', fontsize=9)
    axc.set_title('generic model, ' + r'$\rho = 1$', fontsize=10)

    # --- ring panels: 8 area functions --------------------------------------
    # placed along the ring in cycle order (counter-clockwise from top-left),
    # matching each vowel's angular position in the vowel space:
    # ɨ u o ɔ down the left, a at the bottom, ɛ e i up the right side.
    x = (np.arange(N) + 0.5) / N * L
    slots = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    placed = dict(zip(slots, range(8)))

    for slot, vi in placed.items():
        r, c = slot
        ax = fig.add_axes([0.055 + c * 0.315, 0.715 - r * 0.335, 0.27, 0.27])
        P = mdl.generic_P(thetas_char[vi], rho=1.0, n=N)
        A = mdl.expif(P)
        ax.fill_between(x, 0, A, color='#d0d8f0')
        ax.plot(x, A, 'b-', lw=2)
        ax.plot([0, L], [0, 0], 'b-', lw=2)
        k = indb[vi]
        ax.set_title('[%s]   %s        F1=%d, F2=%d Hz'
                     % (IPA[mdl.VOWELS[vi]], A_PAIRS[vi], f1[k], f2[k]),
                     fontsize=9)
        ax.set_xlim(0, L)
        ax.set_ylim(-0.1, 3.1)
        if r == 2:
            ax.set_xlabel('x (cm)   glottis → lips', fontsize=8)
        else:
            ax.set_xticklabels([])
        if c > 0:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=7)

    savefig(fig, 'Figure1.png')


if __name__ == '__main__':
    main()
