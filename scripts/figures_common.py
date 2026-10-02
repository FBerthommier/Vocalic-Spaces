"""Common helpers for the figure scripts."""

import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
DATA = os.path.join(ROOT, 'data')
FIGDIR = os.path.join(ROOT, 'figures')

sys.path.insert(0, ROOT)


def figdir():
    os.makedirs(FIGDIR, exist_ok=True)
    return FIGDIR


def load_or_compute(npz_name, compute_fn):
    """Load data/<npz_name> if present, else compute and save it."""
    path = os.path.join(DATA, npz_name)
    if os.path.exists(path):
        return dict(np.load(path))
    res = compute_fn()
    np.savez(path, **res)
    return {k: np.asarray(v) for k, v in res.items()}


def load_mc(model, reference=False, quick_n=None):
    """Monte-Carlo C1/C2 data for 'drm' or 'fant'.

    Priority: data/mc_<model>_py.npz (regenerated), then the MATLAB
    reference data/reference/<model>_mc.npz (original draws of the paper).
    """
    py = os.path.join(DATA, 'mc_%s_py.npz' % model)
    ref = os.path.join(DATA, 'reference', '%s_mc.npz' % model)
    if not reference and os.path.exists(py):
        return dict(np.load(py)), 'python'
    if os.path.exists(ref):
        return dict(np.load(ref)), 'matlab-reference'
    raise SystemExit('No Monte-Carlo data for %s: run scripts/monte_carlo.py %s first'
                     % (model, model))


def vowel_space_axes(ax, xlim, ylim):
    """(F2, F1) axes with the F1 axis inverted (Matlab ``axis ij``)."""
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.invert_yaxis()
    ax.set_xlabel('F2 (Hz)')
    ax.set_ylabel('F1 (Hz)')


IPA = {'ɨ': 'ɨ', 'u': 'u', 'o': 'o', 'ɔ': 'ɔ', 'a': 'a', 'ɛ': 'ɛ', 'e': 'e', 'i': 'i'}


def savefig(fig, name):
    path = os.path.join(figdir(), name)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    plt.close(fig)
