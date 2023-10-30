# Atomic-Physics-In-Plasmas-Notes
This repo contains notes on different atomic processes which may occur in a plasma
and toy codes which perform said atomic processes

## Files Descriptions
* Atomic_Collisions_Notes.ipynb - SOS notebook which currently makes use of an R kernel to produce a report on the codes in this repo and accompanying notes
* chargeexchange.py - Python script which performs charge exchange
* unittests.py - Contains the unit tests for the various python scripts in this repo
* style.css - css style file for generated html file from Atomic_Collisions_Notes.ipynb
* Atomic_Collision_Processes_Report.R - Contains R script which is render into a html file
* pre-commit.sh - bash script which must run sucessfully for a commit to be accepted
* run_tests.sh - bash script to run unit tests for repo
* Requirements.txt - Python requirements.txt file used to recreate a reproducible envirnoment

## Requirements
### General programs
* Pandoc (tested with version 3.1.8)
* Python (tested with version 3.11)
* R (tested with version 4.3.1)

### Required R Packages
* knitr (tested with version 1.44)
* bookdown (tested with version 0.35)
* IRdisplay (tested with version 1.1)
* reticulate (tested with version 1.32.0)

### Required Python Packages
* matplotlib (tested with version 3.8.0)
* numpy (tested with version 1.25.2)
* unittest
* scipy (tested with version 1.9.3)
* pygments (tested with version 2.15.1)

### Required Jupyter Kernels
* IR Kernel for R https://github.com/IRkernel/IRkernel (tested with version 1.3.2)
* Sos Kernel https://github.com/vatlab/sos (tested with version 0.24.3)
