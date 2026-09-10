# Green IBIS: Greenland Ice Bed from Ice Surface
Inferring seasonal timescale basal processes from ice surface observations and ice-flow models for Greenland outlet glaciers.

## Note: This repository is actively being built and will undergo many substantial changes.

## Repository Structure
This repository combines MATLAB and Python/Jupyter Notebook scripts to handle observations and models for all parts of the Green IBIS project. 

### data\_tools 
Tools used to obtain datasets for all parts of the project. This folder is divided into distinct tasks that are not already handled by ISSM or other external projects, such as checking for version updates for particular datasets.

### analyses
Discrete analyses that may call on functions in the utils folder.

### runmes
ISSM runme scripts for the project. 

### utils
Utility functions for analyses and pre- and post-processing for models and data.

### plotting
Tools and examples for plotting results.

### notebooks
Jupyter Notebooks.

### tests
Tests for any function that needs to be verfied, particularly in the utils, analyses, and data\_tools folders.
