# PaleOrbitals
Python library to compute and analyze paleo-insolations.

Written in Python2.7

Depends on Python libraries:
  - standard python libraries (math,...)
  - numpy, scipy
  - numba (to accelerate loops through large arrays)

This project is meant to be a simple-to-use Python library with functions
that can compute daily, monthly or annual average insolation values for any
latitude on Earth at any point in time within the range of -50 Myrs to + 20 Myrs
(available through the Laskar2004 orbital parameters dataset: http://www.imcce.fr/Equipes/ASD/insola/earth/earth.html).

It would be nice to have a few methods of computing the same insolation values
to test the method-dependence of output insolations. Right now, the only method
available is based on Laskar2004's fortran program: http://www.imcce.fr/Equipes/ASD/insola/earth/La2004/index.html

Another addition would be functions for the statistical analysis of insolation
output. I'm thinking moslty of like power spectrum time series of insolation as a
function of latitude.But, the main objective of the library is to easily produce
useful insolation values from latitude and time (or orbital parameters directly). The
analyzing and plotting can be done outside this library.

Finally, it would be ideal to have an ipython-notebook that goes though all the calculations in detail.
