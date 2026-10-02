# Matlab supplement (SuppPub1-3) and TLM chain — original files, unmodified

Copies of the supplementary material submitted with the JASA-EL article
("See supplementary material at [URL] for SuppPubN.zip") plus the
transmission-line model used to compute the formants.

## Scripts and their data

| File | Role | Loads |
|---|---|---|
| `simuGEN.m` + `simuGEN.mat` | generic model, Fig. 1 (cycle rho=1, 122 phases) | `.mat`: reference formants f1, f2 and neutral f1b, f2b computed with the compiled TLM |
| `simuDRM.m` + `simuDRM.mat` | DRM, Fig. 3 (C2 cycle + "est"/"SE" estimates) | same principle |
| `simuFANT.m` + `simuFANT.mat` | Fant model, Fig. 4b (C2 cycle) | same principle |
| `TLM.m` | full single-file transmission-line model (P. Badin, ICP): `vtn2frm_ftr_oral` + `nraph_oral` + `aire2spectre*` + `spectrelec` + `peakpick` | — |

**Note on `simuDRM.m`**: this file is byte-identical to the working file
`simuDRM-remplaçant.m` of the author's directory (copy included here under
that name for traceability).  Its first executable line is
`load simuDRM.mat`, so it plots the article figure **without any compiled
TLM**: the formant vectors `f1, f2` and the neutral `f1b, f2b` are read
from the `.mat`, and everything else (coordination function, Fourier
coefficients, "est"/"SE" estimates) is recomputed live.  The actual
`TLM(...)` calls appear as comments, so the same script can be re-run end
to end by uncommenting them with any transmission-line implementation
(the Python equivalent is `../scripts/simu_drm.py`, validated to <1 Hz).

Run in Matlab/Octave from this folder: e.g. `simuDRM` reproduces the
curves of Fig. 3 (`figure(1)`).
