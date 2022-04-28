

This is the code developed for the 2007 thesis
"Numerical Experimentation in Differential Geometry"
by Simon Burton.

https://arrowtheory.com/SimonBurtonThesis.pdf


Pre-requisites
==============

...

Installation
============

...


What's here
===========

Interval Implementations:
-------------------------

This Interval class wraps the mpfi library 
  mpfi.pyx
  mpfi.c
  setup.py
Get libmpfr here: http://www.mpfr.org/mpfr-current/
(NB: needs at least mpfr-2.2.0 with mpfr_get_f patch applied)
Get libmpfi here: http://perso.ens-lyon.fr/nathalie.revol/softwares/

This Interval class solves an optimization problem (useing a more basic Interval class)
  optimize.py

Switch between Interval implementations:
  xinterval.py

Some Tests:
  test_interval.py

Automatic Differentiation
-------------------------

  mjet.py : lazy, multivariate jets

Interface to Taylor ODE Package (Jorba, Zou):
---------------------------------------------
  taylor
  taylor.py
  taylor_tmpl.pyx

Other stuff
-----------

  fmath.py : functions like sqr, sin, etc. that call appropriate methods on arguments

  Makefile
  README.txt


