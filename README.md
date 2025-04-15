# Production of LFR-AS-200 and LFR-AS-30 MPRs to be used in newcleo ORION simulations
This repository contains the necessary scripts to generate the MPR files needed to run ORION simulations of newcleo's LFR-AS-200 and LFR-AS-30 nuclear fuel cycles. The MPR's generated are based on 2024 neutron spectra, which for the LFR-AS-200 is three enrichment zones (inner, middle, outer) and two burnup steps corresponding to the Beginning of Cycle, BoC, and the End of Cycle, EoC. For the LFR-AS-30 the scripts used a single enrichment zone (denoted as Fuel Assembly, FA) and a single burnup step (denoted as single).

## How to download
In order to download this repository the following should be entered into the command line terminal:
```
mkdir MPR_Generation
cd MPR_Generation
git init
git clone git@github.com:sjtaylor-9/LFR-MPR-calculation.git
git pull origin main
```

## How to use
In order to run these scripts a ```conda``` environment with the following packages is required: ```numpy```, ```pandas``` and ```endf```. These can be installed in the command line using either ```pip install endf``` or ```mamba install -c conda-forge endf```, replacing ```endf``` with the respective library name.

This repository is designed so that the user can enter one command into the terminal and all the cross-section data and MPRs will be generated without any need to change the ```python``` logic. This is done by calling a ```shell``` script with ```bash``` in the terminal. When running ```MPR_generation.sh```, three arguments are required. These are:
- The path directory where the output folders/files should be written to,
- The path directory to the local folder containing all of the internal newcleo data, i.e. the neutron flux spectra,
- The path directory to the local folder containing all of the public ENDF/JEFF/TENDL cross-section data.

Here is an example of how to run ```MPR_generation.sh```:
```
bash MPR_generation.sh MPR_Outputs "path/to/folder/with/neutron/fluxes/" "path/to/folder/with/public/cross/section/data/"
```

The ```python``` scripts can also be run individually, noting that an individual set of arguments must be passed in the command line when running. The necessary arguments can be found in the ```parse_arguments()``` function of each script.

**Warnings:**
- Due to the way that the scripts are set-up to read in the neutron flux spectra these must all be stored in the same local folder. The same is true for the ```Lib80x```, ```ace_900```, and ```tendl19c``` folders for the public JEFF, ENDF, and TENDL data sets, respectively.
- If the neutron flux spectra and public cross-section data files are stored in the newcleo OneDrive then the file directories must be given to the ```.sh``` script in "" so that the command line ignores the whitespaces in the directory.
- The python and shell scripts will need to be modified if in the future more burnup steps are to be used.

## Cross-section libraries
The download links for the reaction cross-sections within each of the public libraries are: 
- The data containing the JEFF 3.3 reaction cross-sections can be downloaded from this link: https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/temperatures/ace_900.tar.gz
- The data containing the ENDF reaction cross-sections can be downloaded from this link: https://nucleardata.lanl.gov/lib/Lib80x.zip
- The data containing the TENDL reaction cross-sections can be downloaded from this link: https://tendl.web.psi.ch/tendl_2019/tar_files/tendl19c.tar.bz2

**Warnings:**
- The full unzipped ```Lib80x``` folder containing ENDF data is approximately 20 GB; however, only the ```.ACE``` files within the ```Lib80x``` folder containing the ```.02c``` file extension are needed. This is because these files correlate to a temperature of 900K, which is approximately the temperature of the core. Hence, the rest can be deleted, resulting in a folder of size of 6 GB. Also, some of the file names are incorrectly named, these must unfortunately be manually corrected.

# Authors
These scripts were produced and developed by Micahel Weekes (michael.weekes@newcleo.com) and Sam Taylor (sam.taylor@newcleo.com) **Last modified:** 13/01/2025



