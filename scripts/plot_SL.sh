#!/bin/sh

cp scripts/controlPlotter_SL.tex ${1}/plots_full
cd ${1}/plots_full
pdflatex controlPlotter_SL.tex
cd -
mv ${1}/plots_full/controlPlotter_SL.pdf ${1}.pdf
# pdflatex yields.tex
# cd ../..
