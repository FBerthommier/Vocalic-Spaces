@echo off
rem Build the open-access two-column PDF (MiKTeX / TeXLive)
pdflatex -interaction=nonstopmode paper.tex
bibtex paper
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
