"""Figure 2 -- the coordination function for (a1, a2) (Table I, Eq. 5-7).

The two cosine functions a1(theta) = Psi1 cos(Psi2 - theta) + Omega and
a2(theta) are built from {Omega, Psi1, Psi2} obtained by transforming the
{i, j, k} values; the original points at theta_h in {pi/3, pi, 5 pi/3}
are marked.
"""

import numpy as np
import matplotlib.pyplot as plt

from figures_common import savefig
from vsmodel.models import ijk2psi, ALPHA


def main():
    # {i, j, k} rows for (a1, a2) -- Table I
    ijk = np.array([[1.0, -2.0, 1.0],
                    [1.0, 0.0, -1.0]])
    om, psi1, psi2 = ijk2psi(ijk)
    print('Table I check:  a1: Ω=%g  Ψ1=%g  Ψ2=%g' % (om[0], psi1[0], psi2[0]))
    print('                a2: Ω=%g  Ψ1=%.4f (cos⁻¹π/6=%.4f)  Ψ2=%g'
          % (om[1], psi1[1], 1 / np.cos(np.pi / 6), psi2[1]))

    theta = np.arange(0, 2 * np.pi, 0.01)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    a1 = psi1[0] * np.cos(psi2[0] - theta) + om[0]
    a2 = psi1[1] * np.cos(psi2[1] - theta) + om[1]
    ax.plot(theta, a1, 'b-', lw=1.5, label=r'$a_1(\theta)$')
    ax.plot(theta, a2, 'r-', lw=3, alpha=0.7, label=r'$a_2(\theta)$')

    for ph in (np.pi / 3, np.pi, 5 * np.pi / 3):
        ax.axvline(ph, color='0.5', ls='--', lw=0.8)
    # {i, j, k} points at their phases
    ax.plot([np.pi / 3, np.pi, 5 * np.pi / 3], [1, -2, 1], 'o', color='b', ms=8)
    ax.plot([np.pi / 3, np.pi, 5 * np.pi / 3], [1, 0, -1], 'o', color='r', ms=8)
    ax.annotate('i', (np.pi / 3, 1), xytext=(8, 6), textcoords='offset points', fontsize=12)
    ax.annotate('j', (np.pi, -2), xytext=(8, 6), textcoords='offset points', fontsize=12)
    ax.annotate('k', (5 * np.pi / 3, 1), xytext=(8, 6), textcoords='offset points',
                fontsize=12, color='b')

    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(-3, 3)
    ax.set_xticks([0, np.pi / 3, np.pi / 2, 2 * np.pi / 3, np.pi,
                   4 * np.pi / 3, 3 * np.pi / 2, 5 * np.pi / 3, 2 * np.pi])
    ax.set_xticklabels(['0', 'π/3', 'π/2', '2π/3', 'π', '4π/3', '3π/2', '5π/3', '2π'])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('cosine amplitude')
    ax.legend(loc='lower left')
    ax.grid(alpha=0.25)
    savefig(fig, 'Figure2.png')


if __name__ == '__main__':
    main()
