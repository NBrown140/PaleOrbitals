# PaleOrbitals
Python code to compute and analyze paleo-insolations.

This project is meant to be a simple-to-use Python library with functions
that can compute daily, monthly or annual average insolation values for any
latitude on Earth at any point in time within the range of -50 Myrs to + 20 Myrs
(available through the Laskar2004 dataset).

It would be nice to have a few methods of cumputing the same insolation values
to test the method-dependence of output insolations. Right now, the only method
available is based on Laskar2004's fortran program.
