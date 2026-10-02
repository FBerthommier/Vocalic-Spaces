"""Figure 5 -- 3D view (a1, a2, df1/df2) of conditions C1 and C2:
(a) DRM, (b) Fant model.

Faithful to the submitted figure: the C2 family (coordination function
with random rho, theta) spans two gray surfaces -- the graphs of df1 and
df2 over the (a1, a2) domain -- outlined by the rho = 1 loop (blue: df1,
red: df2); the C1 points (uncorrelated parameters) scatter away from the
surfaces, showing the many-to-one relation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from figures_common import load_or_compute, load_mc, savefig
from simu_drm import drm_cycle
from simu_fant import fant_cycle


def panel(ax, cyc, mc, title, pad=0.12, ngrid=60):
    f1b = float(mc['f1b'])
    f2b = float(mc['f2b'])
    # C2 random draws -> deviations
    df1s = (mc['f1s'] - f1b) / f1b
    df2s = (mc['f2s'] - f2b) / f2b
    # C1 uncorrelated draws
    df1r = (mc['f1r'] - f1b) / f1b
    df2r = (mc['f2r'] - f2b) / f2b
    # rho = 1 loop
    df1c = (cyc['f1'] - f1b) / f1b
    df2c = (cyc['f2'] - f2b) / f2b

    # --- C2 surfaces (graphs of df1 and df2 over the C2 domain) ---------
    x0, x1 = mc['a1s'].min() - pad, mc['a1s'].max() + pad
    y0, y1 = mc['a2s'].min() - pad, mc['a2s'].max() + pad
    gx, gy = np.meshgrid(np.linspace(x0, x1, ngrid),
                         np.linspace(y0, y1, ngrid))
    for zv, gcol in ((df2s, '0.62'), (df1s, '0.78')):
        zg = griddata((mc['a1s'], mc['a2s']), zv, (gx, gy), method='linear')
        ax.plot_surface(gx, gy, zg, color=gcol, shade=True,
                        linewidth=0, antialiased=True, alpha=0.88)

    # --- rho = 1 loop: df1 (blue) and df2 (red) --------------------------
    ax.plot(cyc['a1'], cyc['a2'], df1c, color='b', lw=3, zorder=5)
    ax.plot(cyc['a1'], cyc['a2'], df2c, color='r', lw=3, zorder=5)

    # --- C1 scattered dots -----------------------------------------------
    ax.scatter(mc['a1r'], mc['a2r'], df1r, c='b', s=1.5, alpha=0.45,
               depthshade=False)
    ax.scatter(mc['a1r'], mc['a2r'], df2r, c='r', s=1.5, alpha=0.45,
               depthshade=False)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_zlim(-1, 1)
    ax.set_xlabel(r'$\tilde{a}_1$')
    ax.set_ylabel(r'$\tilde{a}_2$')
    ax.set_zlabel(r'$df_1 / df_2$')
    ax.set_title(title, pad=0)
    ax.view_init(elev=20, azim=-55)
    # legend proxies
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], c='b', lw=3, label=r'$df_1$'),
                       Line2D([], [], c='r', lw=3, label=r'$df_2$')],
              loc='upper right', fontsize=9)


def main(reference=False):
    fig = plt.figure(figsize=(13, 5.6))
    src = ''
    # (a) DRM
    cyc = load_or_compute('drm_cycle_py.npz', lambda: drm_cycle())
    mc, src = load_mc('drm', reference=reference)
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    panel(ax, cyc, mc, '(a) DRM')
    # (b) Fant
    cyc = load_or_compute('fant_cycle_py.npz', lambda: fant_cycle())
    mc, src = load_mc('fant', reference=reference)
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    panel(ax, cyc, mc, '(b) FANT')

    fig.suptitle('C2 coordination surfaces (gray) and rho=1 loop '
                 '(blue df1, red df2) vs C1 dots   [data: %s]' % src,
                 fontsize=10, y=0.99)
    savefig(fig, 'Figure5.png')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', action='store_true',
                    help='use the original MATLAB draws instead of the Python ones')
    main(ap.parse_args().reference)
