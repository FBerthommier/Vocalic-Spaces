"""Convert the MATLAB reference data (.mat) of the article into .npz files.

Run once to (re)build data/reference/*.npz.  The .mat files are read from
the original locations (programs/CompileTLM2024 and programs/); the .npz
copies make the repository self-contained for validation without MATLAB.
"""

import os
import numpy as np
import scipy.io as sio

PROGS = r'D:\HOL\JASA-EL-Modele\newJASA-EL\programs'
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'data', 'reference')
os.makedirs(OUT, exist_ok=True)


def ravel(m, keys):
    return {k: np.ravel(m[k]) for k in keys if k in m}


# cycles computed with the compiled TLM (supplement of the article)
m = sio.loadmat(os.path.join(PROGS, 'CompileTLM2024', 'Supplement', 'simuGEN.mat'))
np.savez(os.path.join(OUT, 'gen_cycle.npz'), **ravel(m, ['f1', 'f2', 'f1b', 'f2b']))

m = sio.loadmat(os.path.join(PROGS, 'CompileTLM2024', 'Supplement', 'simuDRM.mat'))
np.savez(os.path.join(OUT, 'drm_cycle.npz'), **ravel(m, ['f1', 'f2', 'f1b', 'f2b']))

m = sio.loadmat(os.path.join(PROGS, 'CompileTLM2024', 'Supplement', 'simuFANT.mat'))
np.savez(os.path.join(OUT, 'fant_cycle.npz'), **ravel(m, ['f1', 'f2', 'f1b', 'f2b']))

# C1 / C2 Monte-Carlo data of Figures 3, 5, 6, 7
m = sio.loadmat(os.path.join(PROGS, 'simurandDRM4expcos.mat'))
np.savez(os.path.join(OUT, 'drm_mc.npz'), **ravel(
    m, ['a1', 'a2', 'a10', 'a20', 'f1', 'f2', 'f1est', 'f2est', 'Pval',
        'indthetab', 'f1b', 'f2b',
        'f1r', 'f2r', 'a1r', 'a2r', 'a10r', 'a20r',
        'f1s', 'f2s', 'a1s', 'a2s', 'a10s', 'a20s', 'rhos', 'thetas']))

m = sio.loadmat(os.path.join(PROGS, 'simurandFANT2.mat'))
np.savez(os.path.join(OUT, 'fant_mc.npz'), **ravel(
    m, ['a1', 'a2', 'a10', 'a20', 'f1', 'f2', 'Pval', 'indthetab', 'f1b', 'f2b',
        'f1r', 'f2r', 'a1r', 'a2r',
        'f1s', 'f2s', 'a1s', 'a2s', 'rhos', 'thetas']))

# glottal excitation used by the LPC resynthesis (ASCII 'excit.mat')
exc = np.loadtxt(os.path.join(PROGS, 'CompileTLM2024', 'excit.mat'))
np.save(os.path.join(OUT, 'excit.npy'), exc)

print('reference data written to', OUT)
for f in sorted(os.listdir(OUT)):
    print('  ', f)
