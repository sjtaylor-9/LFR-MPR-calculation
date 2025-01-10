# Production of LFR-AS-200 and LFR-AS-30 MPRs to be used in newcleo ORION simulations
This repository contains the necessary scripts to generate the MPR files needed to run ORION simulations of newcleo's LFR-AS-200 and LFR-AS-30 nuclear fuel cycles. The MPR's generated are based on 2024 neutron spectra, which for the LFR-AS-200 is 3 enrichment zones (inner, middle, outer) and 2 burnup steps corresponding to the Beginning of Cycle, BoC, and the End of Cycle, EoC. For the LFR-AS-30 the scripts used a single enrichment zone (denoted as Fuel Assembly, FA) and a single burnup step (denoted as single).

## How to use
In order to download this repository the following should be entered into the command line terminal:
```
mkdir MPR_Generation
cd MPR_Generation
git init
git clone git@github.com:sjtaylor-9/LHCb-D0-asymmetry-2.git
git pull origin main
```

In order to run these scripts a conda environment with the following packages is required: ```numpy```, ```pandas``` and ```endf```. These can be installed in the command line using either ```pip install $x$``` or if ```mamba install -c conda-forge $x$```, where $x$ is the name of the python library.

# Authors
These scripts were produced and developed by Micahel Weekes (michael.weekes@newcleo.com) and Sam Taylor (sam.taylor@newcleo.com) **Last modified:** 10/01/2025



