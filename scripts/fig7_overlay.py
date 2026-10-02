"""Figure 7 -- DRM C2 overlaid on the human vowel space ("Adult males",
after Boe et al. 2019, Fig. 4, L = 17.5 cm), reproduced from the
submitted figure.

The human-space envelope, the 7-parameter Monte-Carlo cloud and the
vowel ellipses are DIGITIZED from the submitted Figure7.jpg (calibration
by the axis gridlines; see tools/digitize_fig7.py).  The DRM C2 loop is
computed by this repository (red circles at the eight characteristic
phases).  Axes in kHz with the F2 axis reversed, as in the original.
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from figures_common import DATA, load_or_compute, savefig
from simu_drm import drm_cycle

HERE = os.path.dirname(os.path.abspath(__file__))
DIG = os.path.join(DATA, 'fig7_digitized.npz')   # digitized from the
# submitted Figure7.jpg (Boe et al. 2019 overlay); build with
# tools/digitize_fig7.py (calibration by the axis gridlines)

# vowel ellipses digitized / read on the submitted figure:
# (F2, F1, semi-F2, semi-F1) in kHz, color, label
ELLIPSES = [
    (2.29, 0.285, 0.185, 0.062, 'red',       'i'),
    (2.09, 0.415, 0.100, 0.075, 'red',       'I'),
    (1.99, 0.545, 0.085, 0.065, 'red',       r'$\varepsilon$'),
    (0.78, 0.305, 0.130, 0.075, 'green',     'u'),
    (0.90, 0.425, 0.095, 0.065, 'green',     'U'),
    (0.72, 0.555, 0.130, 0.075, 'c',         'o'),
    (1.07, 0.625, 0.100, 0.065, 'navy',      r'$\Lambda$'),
    (0.83, 0.705, 0.060, 0.050, 'navy',      'O'),
    (1.42, 0.745, 0.070, 0.050, 'navy',      r'$\alpha$'),
]
ECOL = {'red': 'red', 'green': 'limegreen', 'c': 'c', 'navy': 'navy'}


def main(reference=False):
    fig, ax = plt.subplots(figsize=(8.6, 7.6))

    # ---- digitized human space -----------------------------------------
    if os.path.exists(DIG):
        d = np.load(DIG, allow_pickle=True)
        env_f2, env_f1 = d['env_f2'], d['env_f1']
        ax.fill(env_f2, env_f1, color='0.70', zorder=1)
        ax.plot(env_f2, env_f1, 'k-', lw=1.6, zorder=3)
        ax.plot(d['cloud_f2'], d['cloud_f1'], '.', color='b', ms=2.2,
                alpha=0.55, zorder=2, rasterized=True)
    else:
        raise SystemExit('missing %s : run tools/digitize_fig7.py' % DIG)

    # ---- vowel ellipses --------------------------------------------------
    for f2, f1, w2, h1, col, lab in ELLIPSES:
        ax.add_patch(Ellipse((f2, f1), 2 * w2, 2 * h1, fill=False,
                             edgecolor=ECOL[col], lw=1.6, zorder=4))
        ax.text(f2, f1, lab, fontsize=15, ha='center', va='center',
                color=ECOL[col], zorder=5)

    # ---- DRM C2 (computed here): red circles at the 8 vowel phases ------
    d = load_or_compute('drm_cycle_py.npz', lambda: drm_cycle())
    indb = np.array([1, 21, 31, 41, 61, 81, 91, 101]) - 1
    ax.plot(d['f2'][indb] / 1000.0, d['f1'][indb] / 1000.0, 'o', ms=5,
            mfc='none', mec='red', mew=1.2, ls='none', zorder=6)

    # ---- axes: kHz, F2 reversed, F1 inverted (0.2 kHz on top) -----------
    ax.set_xlim(2.7, 0.2)
    ax.set_ylim(1.03, 0.12)
    ax.set_xticks([2.5, 2.0, 1.5, 1.0, 0.5])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel('F2 (kHz)')
    ax.set_ylabel('F1 (kHz)')
    ax.set_title('Adult males')
    ax.grid(color='0.85', lw=0.6, zorder=0)
    savefig(fig, 'Figure7.png')


if __name__ == '__main__':
    main()
