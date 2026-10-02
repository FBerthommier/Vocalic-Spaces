"""Multimedia files Mm. 1-6 of the article, regenerated in Python.

* Mm. 1/4/6 (.wav): the 8 vowel sounds at the characteristic phases for
  the generic, DRM and Fant models (LPC resynthesis from the TLM spectra,
  0.5 s per vowel + 75 ms silence, 20 kHz).
* Mm. 2/3/5 (.gif, .mp4 if ffmpeg is available): the area function
  A(x, theta) from 0 to 2 pi with the running point in the vowel space.

The originals (submitted to JASA-EL) are kept in multimedia/ for
comparison.
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

from vsmodel import tlm, models as mdl
from vsmodel.synthesis import voysynth, wavwrite

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
MM = os.path.join(ROOT, 'multimedia', 'generated')
EXCIT = np.load(os.path.join(ROOT, 'data', 'reference', 'excit.npy'))


def synth_8_wavs(name, specs):
    sig = []
    for spec in specs:
        sig.append(voysynth(spec, excit=EXCIT))
        sig.append(np.zeros(1500))
    out = os.path.join(MM, name)
    wavwrite(np.concatenate(sig), 20000, out)
    print('saved', out)


def make_movie(name, areas, f1s, f2s, L, n, step_frame=1):
    """Movie: area function (top) + running point in the vowel space (bottom)."""
    x = (np.arange(n) + 0.5) / n * L
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 6.4))
    ax1.set_xlim(0, L)
    ax1.set_ylim(0, max(np.max(a) for a in areas) * 1.05)
    ax1.set_xlabel('x (cm)   glottis → lips')
    ax1.set_ylabel('A (cm²)')
    line1, = ax1.plot([], [], 'b-', lw=2)
    ax1.fill_between(x, 0, areas[0], color='#d0d8f0')
    ax2.plot(f2s, f1s, 'b-', lw=0.8)
    ax2.set_xlabel('F2 (Hz)')
    ax2.set_ylabel('F1 (Hz)')
    ax2.invert_yaxis()          # article convention (Matlab axis ij):
                                # low F1 at top -- the triangle is not flipped
    point2, = ax2.plot([], [], 'ro', ms=8)

    def update(i):
        line1.set_data(x, areas[i])
        for coll in list(ax1.collections):
            coll.remove()
        ax1.fill_between(x, 0, areas[i], color='#d0d8f0')
        point2.set_data([f2s[i]], [f1s[i]])
        return line1, point2

    anim = FuncAnimation(fig, update, frames=range(0, len(areas), step_frame),
                         blit=False)
    gif = os.path.join(MM, name + '.gif')
    anim.save(gif, writer=PillowWriter(fps=12))
    print('saved', gif)
    try:
        mp4 = os.path.join(MM, name + '.mp4')
        anim.save(mp4, writer=FFMpegWriter(fps=12))
        print('saved', mp4)
    except Exception as e:
        print('(mp4 skipped: %s)' % e)
    plt.close(fig)


# the 8 characteristic phases (same for the three models: 0, pi/3, pi/2,
# 2pi/3, pi, 4pi/3, 3pi/2, 5pi/3)
THETAS8 = np.pi * np.array([0, 1 / 3, 1 / 2, 2 / 3, 1, 4 / 3, 3 / 2, 5 / 3])


def gen_all(quick=False):
    os.makedirs(MM, exist_ok=True)
    pastheta = np.pi / 15 if quick else np.pi / 60
    step_frame = 2 if quick else 1
    ntheta = int(round(2 * np.pi / pastheta)) + 2
    theta = np.linspace(0, pastheta * (ntheta - 1), ntheta)

    # ---------------- generic (Mm. 1, Mm. 2) ------------------------------
    t0 = time.time()
    areas, f1s, f2s, specs8 = [], [], [], []
    for th in theta:
        P = mdl.generic_P(th, 1.0, n=100)
        f1, f2, f3, spec = tlm.TLM3spec(mdl.generic_area(P), 1)
        areas.append(mdl.expif(P)); f1s.append(f1); f2s.append(f2)
    for th in THETAS8:                        # the 8 characteristic vowels
        P = mdl.generic_P(th, 1.0, n=100)
        specs8.append(tlm.TLM3spec(mdl.generic_area(P), 1)[3])
    synth_8_wavs('MM1.wav', specs8)
    make_movie('MM2', areas, f1s, f2s, L=17.5, n=100, step_frame=step_frame)
    print('generic done in %.0f s' % (time.time() - t0))

    # ---------------- DRM (Mm. 3, Mm. 4) -----------------------------------
    Om, Ps1, Ps2 = mdl.drm_setup()
    areas, f1s, f2s, specs8 = [], [], [], []
    for th in theta:
        Pval = Om + Ps1 * np.cos(Ps2 - th)
        area = mdl.drm_area(Pval, n=120)
        f1, f2, f3, spec = tlm.TLM3spec(area, 1)
        areas.append(area[:, 1]); f1s.append(f1); f2s.append(f2)
    for th in THETAS8:
        Pval = Om + Ps1 * np.cos(Ps2 - th)
        specs8.append(tlm.TLM3spec(mdl.drm_area(Pval, n=120), 1)[3])
    synth_8_wavs('MM4.wav', specs8)
    make_movie('MM3', areas, f1s, f2s, L=17.5, n=120, step_frame=step_frame)

    # ---------------- Fant (Mm. 5, Mm. 6) ----------------------------------
    Om, Ps1, Ps2, Lc, l, Lip, A = mdl.fant_setup(200)
    areas, f1s, f2s, specs8 = [], [], [], []
    for th in theta:
        Pval = Om + Ps1 * np.cos(Ps2 - th)
        area = mdl.fant_area(Pval, Lc, l, Lip, A, n=200)
        f1, f2, f3, spec = tlm.TLM3spec(area, 0)
        areas.append(area[:, 1]); f1s.append(f1); f2s.append(f2)
    for th in THETAS8:
        Pval = Om + Ps1 * np.cos(Ps2 - th)
        specs8.append(tlm.TLM3spec(mdl.fant_area(Pval, Lc, l, Lip, A, n=200), 0)[3])
    synth_8_wavs('MM6.wav', specs8)
    make_movie('MM5', areas, f1s, f2s, L=17.5, n=200, step_frame=step_frame)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true',
                    help='fewer frames (fast preview, lower quality)')
    args = ap.parse_args()
    gen_all(args.quick)
