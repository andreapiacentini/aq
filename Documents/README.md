# AQ
Air Quality 4DEnVar suite for the SEEDS project based on the JEDI Data Assimilation environment

## Context and Credits
Talk about the SEEDS projects, its fundings, but also the JEDI project at JCSDA and the ECMWF work on Atlas and on the qg test case used as template

## Position of the problem
Quick overview on the theory

### Air Quality modelling
Something about the models, what in the state, the internal coherency, the assumed sources of uncertainty

### Observations for Air Quality
At least the ground stations we are using now

### Ensemble Data Assimilation
Not a full course, but a minimum of wording that could help understanding the following. In particular the localization.

## OOPS and SABER
The usual few lines of presentation.  
How they help us and the idea of a "coding contract" for the interfaces.

### Main classes
If we want to go to this detail

### Yaml input files
Just the principle and a link to an explanation of the yaml syntax. No need to detail them one by one.

## Actual implementation
What in our specific 4DEnVar

### Model ensemble
Our geometry on a rectangular regular limited area grid.  
The NetCDF storage with its loader and writer.

### Observations
The HDF5 based storage derived from DAIMON with its loader and writer

### The algorithm
Features of our 4DEnVar.  
Multivariate localization.  
NO2 log transformation.  
Quality control filters on the initial misfit.  
Bias correction.

## Known issues and future developments
Let's see if we have to mention here the not so canonical interpolation, the need of optimization of the localization in view of the increase in size and resolution of the domain, the advection of the localization and what else.
