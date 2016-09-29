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


**Reference:** Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015. Accurate polynomial expressions for the 
density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43.


## Fortran code

The file `polyTEOS10_bsq.F90` includes fortran procedures to compute in-situ density, thermal expansion and haline contraction coefficients.
It also provides an approximation to convert conservative temperature into potential temperature, useful in ocean models
to compute air-sea fluxes.

We propose only the **polyTEOS10_bsq** approximation in Fortran, as this is the one that is generally most well
suited for most ocean models. Note that the ocean model [NEMO](http://www.nemo-ocean.eu/) uses this equation of state
as the standard since version 3.6.

A fortran version of **polyTEOS10_75t** is also available as part of the standard [TEOS-10 package](http://www.teos-10.org/software.htm).


## Python code

The four versions of polyTEOS10 are coded in the python script `polyTEOS10.py`.

## Matlab code

Each version of polyTEOS10 is also distributed for Matlab:
* `polyTEOS10_bsq.m`
* `polyTEOS10_55t.m`
* `polyTEOS10_75t.m`
* `polyTEOS10_stif.m`

## Contributors

Fabien Roquet (Department of Meteorology at Stockholm University). 
More information on the author on his [personal webpage](http://www.su.se/profiles/froqu)


## License

This software is distributed under the terms of the GNU GENERAL PUBLIC LICENSE v3. The GNU General Public License is a free, copyleft license for software and other kinds of works.



