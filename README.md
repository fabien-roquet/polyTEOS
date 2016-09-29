# polyTEOS
Polynomial versions of TEOS-10 seawater equation of state (Fortran, Matlab and Python)

## Motivation

A new set of approximations to the standard TEOS-10 equation of state has recently been proposed. 
These follow a polynomial form, making it computationally efficient for use in numerical ocean models. 
Two versions are provided, the first **polyTEOS10\_bsq** being a fit of density for Boussinesq ocean models, 
and the second **polyTEOS10\_55t** fitting specific volume which is more suitable for compressible models.

Both versions are given as the sum of a vertical reference profile (6th-order polynomial) 
and an anomaly (52-term polynomial, cubic in pressure), with relative errors of 0.1% on the 
thermal expansion coefficients.

A 75-term polynomial expression **polyTEOS10_75t** is also presented for computing 
specific volume, with a better accuracy than the existing TEOS-10 48-term rational approximation, 
especially regarding the sound speed. This version is now used (since 2015) as the standard 
approximation of the [TEOS-10 equation of state](http://www.teos-10.org/).

Finally, a stiffened version **polyTEOS10\_stif** of the equation of state is also proposed.

More details on these approximations of the TEOS-10 equation of state are available in Roquet et al. (2016).


Reference: Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015. Accurate polynomial expressions for the 
density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43.

Calculates in-situ density from Absolute Salinity, Conservative
Temperature and pressure, using the computationally-efficient 55-term
polynomial expression for density (Roquet et al., 2014).


Note that the 55-term equation has been fitted in a restricted range of
parameter space, and is most accurate inside the "oceanographic funnel"
described in McDougall et al. (2003).  The GSW library function
"gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
some of one's data lies outside this "funnel".




