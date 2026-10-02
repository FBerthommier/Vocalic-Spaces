# Open-access two-column edition of the article

This folder holds the **LaTeX source, the compiled `paper.pdf` (11
pages) and the arXiv submission package** of *A mathematical model of
the vowel space*, F. Berthommier, in a two-column layout equivalent to
the JASA-EL "reprint" form, built exclusively with free components
(standard `article` class, `natbib`/`plainnat`, `TIPA`, `newtx`,
`empheq`, `caption`, `hyperref`).

The canonical, citable record is **arXiv:2111.00868v3**
(<https://arxiv.org/abs/2111.00868v3>; note: v2 is a dead version).
A didactic introduction is provided in
[`../docs/tutorial.pdf`](../docs/tutorial.pdf).

## Overleaf / arXiv

* **arXiv (v3)**: upload `arxiv_submission.zip` as-is.  It contains the
  flat set `paper.tex`, `Paperbib.bib`, the **pre-compiled `paper.bbl`**
  (arXiv does not run BibTeX) and the seven figures.
* **Overleaf**: create a new project from the same files; compile with
  pdfLaTeX.
* Build locally: `build.bat` (pdfLaTeX + BibTeX + 2 passes) → produces
  `paper.pdf`.

The repository URL is set in `paper.tex`:

```latex
\newcommand{\repobase}{https://github.com/FBerthommier/Vocalic-Spaces}
```

## Cosmetic changes applied to the submission source

Body text, equations, tables, captions and bibliography are carried over
**verbatim** from `programs/Resoumission/Resubmitted/Paper.tex`.  The
changes are limited to:

1. **template swap** — `JASA.cls` → free `article[twocolumn]` stack
   (header, title block, `\docsection`→`\section*`,
   `\multimedia`→linked list, `\linenomath`/`\reprintcolumnwidth` no-ops,
   `jasaauthyear2.bst`→`plainnat`);
2. **figures centered** (`\centering` in the seven figure environments —
   verified programmatically at 0.00 pt offset);
3. **running title + arXiv reference on every page**, including the
   title page (fancyhdr `plain` redefinition);
4. caption justification via the `caption` package, tables typeset
   `\small` with `footnotesize` captions;
5. footnotes "[URL will be inserted by AIP] SuppPubN.zip" → links to the
   companion repository; submission-only end matter (duplicated
   "Figure Captions" list) dropped;
6. Table III column spec `{ccccc}`→`{cccc}` (the original declared 5
   columns for 4 used; no content change).

### Layout fixes (v2)

* **p. 6 blank removed**: `\dbltextfloatsep` 20→10 pt and reduced
  display skips around Eq. 8 pull the equation into the bottom of the
  right column;
* **p. 9 mid-column blank removed**: `\raggedbottom`;
* **single "References" heading** (`plainnat` emits its own);
* **Tables page**: `\onecolumn` break, deterministic 2×2 grid (two
  minipages, captions via `\captionof`) — tables I+II left, III+IV right;
* **Table II enlarged**: `\arraystretch{1.7}` so the four
  `cos⁻¹(atan(...))` formulas never touch the rules;
* Figure 3 slightly reduced (12→11.7 cm) to balance its page.

## Strict verification

`../scripts/verify_paper.py` compares this edition with the original
submission source and checks that nothing else changed:

```
python scripts/verify_paper.py [path/to/original/Paper.tex]
```

Current result: body text 0 unexpected differences; equations 13/13;
table rows 17/17; figure captions 7/7; Mm. descriptions 6/6; phonetic
strings complete.  The standalone `arxiv_submission.zip` was compiled in
isolation (no BibTeX, as on arXiv): 11 pages, 0 errors, 0 unresolved
citations.

## Files

- `paper.tex` — source (free template)
- `Paperbib.bib` — bibliography (unchanged from the submission)
- `paper.bbl` — pre-compiled bibliography for arXiv
- `Figure1.jpg` … `Figure7.jpg` — figures (unchanged, author's files)
- `arxiv_submission.zip` — flat Overleaf/arXiv package
- `build.bat` — local build script

## License

The **whole repository — including this article edition (LaTeX source,
figures, multimedia) — is distributed under the MIT License**
(top-level `LICENSE`), with authorship credits to Frédéric Berthommier
(article, simulations, package), **Pierre Badin** (transmission-line
model of the vocal tract, original MATLAB codes) and **Laurent Girin**
(glottal source and LPC resynthesis, original MATLAB codes).  The
compiled article states this in its `\thanks` note; the canonical
citable record is arXiv:2111.00868.
