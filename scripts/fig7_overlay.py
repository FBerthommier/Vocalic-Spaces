"""Figure 7 -- overlay of the DRM C2 data with the human vowel space
(ellipses at 2 sigma, Peterson & Barney 1952 male speakers).

The article overlays it with the 7-parameter Monte-Carlo cloud of
Boe et al. (2019, Fig. 4, L = 17.5 cm, Science Advances 5(12), eaaw3916,
open access).  That cloud is not redistributable; if you digitize it into
``data/boe2019_fig4.csv`` (two columns f1,f2 in Hz), it is drawn in blue
and the figure matches the published one.  Without the file, the DRM C2
domain and the human vowel space are still reproduced.

Only the [i] area of the human vowel space is not covered by the DRM C2
domain (2-parameter model).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from figures_common import DATA, load_or_compute, load_mc, savefig
from simu_drm import drm_cycle

# Peterson & Barney (1952), 76 male speakers: mean F1/F2 (Hz) and SDs
PB1952 = {
    'i':  (270, 2290, 58, 204),
    'ɪ':  (390, 1990, 71, 306),
    'ɛ':  (530, 1840, 87, 394),
    'æ':  (660, 1720, 95, 349),
    'ɑ':  (730, 1090, 103, 328),
    'ɔ':  (570, 840, 85, 313),
    'ʊ':  (440, 1020, 90, 438),
    'u':  (300, 870, 54, 321),
    'ʌ':  (640, 1190, 105, 383),
    'ɝ':  (490, 1350, 88, 341),
}


def draw_ellipses(ax, k=2.0):
    for lab, (f1, f2, s1, s2) in PB1952.items():
        ax.add_patch(Ellipse((f2, f1), 2 * k * s2, 2 * k * s1,
                             facecolor='0.85', edgecolor='0.55', lw=0.8, zorder=1))
        ax.annotate(lab, (f2, f1), fontsize=11, ha='center', va='center',
                    color='0.35', zorder=2)


def load_boe2019():
    """Optional digitized cloud of Boe et al. (2019), Fig. 4 (L = 17.5 cm)."""
    path = os.path.join(DATA, 'boe2019_fig4.csv')
    if os.path.exists(path):
        d = np.loadtxt(path, delimiter=',', skiprows=1)
        return d[:, 0], d[:, 1]
    return None, None


def main(reference=False):
    d = load_or_compute('drm_cycle_py.npz', lambda: drm_cycle())
    mc, src = load_mc('drm', reference=reference)

    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    draw_ellipses(ax)

    # optional: digitized cloud of Boe et al. (2019) 7-parameter Monte-Carlo
    f1_boe, f2_boe = load_boe2019()
    note = ''
    if f1_boe is not None:
        ax.plot(f2_boe, f1_boe, '.', color='cornflowerblue', ms=1.2, alpha=0.4,
                zorder=3, label='Boë et al. (2019), Fig. 4, L=17.5 cm')
    else:
        note = '  (Boë et al. 2019 cloud: see docstring)'

    # DRM C2 domain: random (rho, theta) configurations + rho=1 cycle
    ax.plot(mc['f2s'], mc['f1s'], '.', color='darkgreen', ms=1.2, alpha=0.5,
            zorder=4, label='DRM C2 (rho, theta random)')
    ax.plot(d['f2'], d['f1'], 'r-', lw=1.8, zorder=5, label='DRM C2, rho=1')
    ax.plot(d['f2'][[20, 60, 100]], d['f1'][[20, 60, 100]], 'o', color='r',
            ms=5, zorder=6)
    ax.annotate('[i] not covered', (2290, 270), xytext=(1500, 250),
                fontsize=10, arrowprops=dict(arrowstyle='->', color='0.3'))

    ax.set_xlim(500, 2800)
    ax.set_ylim(100, 900)
    ax.invert_yaxis()
    ax.set_xlabel('F2 (Hz)')
    ax.set_ylabel('F1 (Hz)')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title('DRM C2 vs human vowel space [C2 data: %s]%s' % (src, note),
                 fontsize=10)
    savefig(fig, 'Figure7.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB draws instead of the Python ones')
    main(ap.parse_args().reference)
