"""Digitize the submitted Figure7.jpg (Boe et al. 2019 overlay).

Builds ``../data/fig7_digitized.npz`` used by :mod:`fig7_overlay`:
  * human-space envelope (black outline of the gray region),
  * 7-parameter Monte-Carlo cloud (blue pixels, subsampled),
  * axis calibration by the gridlines (F2 reversed: 328 px -> 2.5 kHz,
    168 px / 0.5 kHz; F1: 222 px -> 0.2 kHz, 181 px / 0.2 kHz).
The vowel ellipses used by fig7_overlay.py are read visually on the
original and hard-coded there.  Run once from this folder:

    python digitize_fig7.py
"""

import os

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, '..', 'article', 'Figure7.jpg')
OUT = os.path.join(HERE, '..', 'data', 'fig7_digitized.npz')

img = np.asarray(Image.open(IMG).convert('RGB')).astype(int)
H, W = img.shape[:2]
R, G, B = img[..., 0], img[..., 1], img[..., 2]


def F2(x):
    return 2.5 - (x - 328) * 0.5 / 168.0


def F1(y):
    return 0.2 + (y - 222) * 0.2 / 181.0


BX0, BX1, BY0, BY1 = 265, 1099, 137, 941
out = {}

# ---- envelope: largest black component inside the box ------------------
black = (R < 90) & (G < 90) & (B < 90)
inner = black.copy()
inner[:BY0, :] = False
inner[BY1:, :] = False
inner[:, :BX0] = False
inner[:, BX1:] = False
lab, n = ndimage.label(inner)
sizes = ndimage.sum(inner, lab, range(1, n + 1))
ys, xs = np.where(lab == (np.argmax(sizes) + 1))
cx, cy = xs.mean(), ys.mean()
order = np.argsort(np.arctan2(ys - cy, xs - cx))
out['env_f2'] = F2(xs[order][::8])
out['env_f1'] = F1(ys[order][::8])

# ---- blue cloud: subsampled pixels --------------------------------------
blue = (B > 140) & (B > R + 40) & (G < B - 20)
blue[:BY0, :] = False
blue[BY1:, :] = False
blue[:, :BX0] = False
blue[:, BX1:] = False
ys, xs = np.where(blue)
sel = np.random.RandomState(0).choice(len(xs),
                                      size=min(9000, len(xs)), replace=False)
out['cloud_f2'] = F2(xs[sel])
out['cloud_f1'] = F1(ys[sel])

np.savez(OUT, **out)
print('written', OUT, '| envelope %d pts, cloud %d pts'
      % (len(out['env_f2']), len(out['cloud_f2'])))
