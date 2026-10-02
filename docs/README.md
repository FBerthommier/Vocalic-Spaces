# Documentation

## The Vowel Space Model — A Didactic Introduction

`tutorial.tex` / `tutorial.pdf` (7 pages, US English) — a step-by-step
didactic development of the theory of the research article:

* the many-to-one articulatory–acoustic problem;
* the Schroeder–Ehrenfest perturbation relation and the odd cosine
  Fourier coefficients;
* the three-phase mixing function (Gibbs triangle → cycle) and the seed
  vectors {i, j, k};
* the generic model, its soft rectifier and the eight cardinal vowels
  (with the formant table computed by this repository's validated code);
* the coordination function {Ω, Ψ1, Ψ2} and its 4-tube applications
  (DRM, Fant), conditions C1/C2 and the bijection result;
* the quadratic bias of the Schroeder–Ehrenfest relation;
* a guided tour of the repository scripts, exercises, and pointers.

Build: `pdflatex tutorial` (twice).  Uses the same free toolchain as
`../paper/` (article + natbib + tipa + newtx).

## Pointers to the research article

* **arXiv (canonical)**: <https://arxiv.org/abs/2111.00868> —
  *A mathematical model of the vowel space*, F. Berthommier.
* The article itself (PDF) is intentionally **not** hosted in this
  repository; cite it through the arXiv record (see `../CITATION.cff`).
* Open-access LaTeX source of the article: `../paper/` (compile or use
  `../paper/arxiv_submission.zip`).

## License

This documentation and the whole repository are distributed under the
**MIT License** (see `../LICENSE`): Frédéric Berthommier (article,
simulations, package), Pierre Badin (transmission-line model, original
MATLAB), Laurent Girin (glottal source and LPC resynthesis, original
MATLAB).
