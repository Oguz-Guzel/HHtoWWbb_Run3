#!/bin/sh

cp scripts/controlPlotter_DL.tex ${1}/plots
cd ${1}/plots
pdflatex controlPlotter_DL.tex
cd -
mv ${1}/plots/controlPlotter_DL.pdf ${1}.pdf
# pdflatex yields.tex
# cd ../..
