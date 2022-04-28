all: pydx/scalar/mpfi.so

pydx/scalar/mpfi.so: pydx/scalar/mpfi.c
	python setup.py build_ext -i

pydx/scalar/mpfi.c: pydx/scalar/mpfi.pyx
	pyrexc pydx/scalar/mpfi.pyx

