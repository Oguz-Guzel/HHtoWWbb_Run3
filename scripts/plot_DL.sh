#!/bin/sh

cp scripts/controlPlotter_DL.tex ${1}/plots_full
cd ${1}/plots_full
pdflatex controlPlotter_DL.tex
cd -
mv ${1}/plots_full/controlPlotter_DL.pdf ${1}.pdf
# pdflatex yields.tex
# cd ../..
