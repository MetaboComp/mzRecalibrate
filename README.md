# MzRecalibrate

## To do
- Uncalibrated data?

**Recalibration of m/z axis in mzML objects (e.g. from LC-HRMS)**  
Associate Professor Carl Brunius  <carl.brunius@chalmers.se>  
Department of Life Sciences
Chalmers University of Technology www.chalmers.se

## General description
The package contains functions to identify "dragging" m/z traces
suitable for mass axis (m/z) calibration, to visualize these traces for manual curation
and to perform actual mass axis recalibration on a scan by scan basis.

## Installation
- Please ensure that R is installed (https://www.r-project.org/)
- I also recommend installing RStudio for a smoother R experience (https://rstudio.com/) 

mzRecalibrate depends on some packages from BioConductor
```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("mzR", quietly = TRUE)) BiocManager::install("mzR")
if (!require("MSnbase", quietly = TRUE)) BiocManager::install("MSnbase")
if (!require("xcms", quietly = TRUE)) BiocManager::install("xcms")
```

To install the mzRecalibrate package you also need the `remotes` R package
```
if (!require("remotes", quietly = TRUE)) install.packages('remotes')
remotes::install_gitlab('CarlBrunius/mzRecalibrate')
```

## Version history
version | date | comment
:------ | :--- | :------
0.1.00  | 2023-05-25 | Hello world. Transferred function files from early development.
