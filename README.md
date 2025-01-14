![](https://img.shields.io/badge/CMS-Run3-blue)

# HH->WWbb Run-3 analysis
----WORK IN PROGRESS----

This repository uses the **bamboo analysis framework**, you can install it via the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html#fresh-install

and install CMSJMECalculators with correctionlib
```sh
git clone https://gitlab.cern.ch/cp3-cms/CMSJMECalculators.git
pip install ./CMSJMECalculators
pip install correctionlib
```

Then clone this repository in the parent directory containing the bamboo installation:

```sh
git clone ssh://git@gitlab.cern.ch:7999/aguzel/HHtoWWbb_Run3.git && cd HHtoWWbb_Run3
```

Execute these each time you start from a clean shell on lxplus or any other machine with an cvmfs:
```sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
source (path to your bamboo installation)/bamboovenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:${PWD}/python/"
```

and the followings before submitting jobs to the batch system (HTCondor, Slurm, Dask and Spark are supported):

```sh
voms-proxy-init --voms cms -rfc --valid 192:00 
export X509_USER_PROXY=$(voms-proxy-info -path)
```
if you encounter problems with accessing files when using batch, the following lines may solve your problem

```sh
voms-proxy-init --voms cms -rfc --valid 192:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Then cutflow study of the analysis is executed via the following command line using batch (you can pass `--maxFiles 1` to use only 1 file from each sample for a quick test):

```sh
bambooRun -m python/cutflowAnalysis.py config/<2022 or 2023>_v12.yml -o ./outputDir/ --envConfig config/cern.ini  --distributed driver
```
Instead of passing `--envConfig config/cern.ini` everytime, you can copy the content of that file to `~/.config/bamboorc`.

using the skims in the root files that re in the results directory of bamboo output `./outputDir/`, you can perform machine learning applications. Once you get the output from your machine learning (preferably in onnx format), you can evaluate that on the datasets by

```sh
bambooRun -m python/mvaEvaluator.py config/<2022 or 2023>_v12.yml -o ./outputDir/ --envConfig config/cern.ini  --distributed driver
```
where you get the ML score applied on the analysis.

To produce sync skim
```sh
bambooRun -m python/syncSkimmer.py config/2022_v12_sync.yml -o output/syncTest --sync
```