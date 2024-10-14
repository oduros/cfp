# Coefficients of fractional parentage

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10700966.svg)](https://doi.org/10.5281/zenodo.10700966)

We made available an open source database of coefficients of fractional parentage (cfp) compiled from Dobromir D. Velkov PhD thesis. The database was constructed by a Python 3 script designed to extract consistently the listed tables from the source PDF. 

## Database
The database is accessible here: https://doi.org/10.5281/zenodo.10700966
It follows the listing conventions introduced by Racah and Nielson & Koster. 

## Script
Here are provided scripts to use this database, one for computing matrix elements of Racah unit tensor operator, one to apply its usage for the development of the crystal electric field. In the latter, an example with Er3+ in a tetragonal environment is provided.
Requirements :
 * Python 3
 * Numpy
 * Sympy
 * h5py
