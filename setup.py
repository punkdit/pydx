from distutils.core import setup, Extension
import distutils.sysconfig

#vars = distutils.sysconfig.get_config_vars()
#vars['OPT'] = ' -w -g ' #-O3

ext = Extension(
    "pydx.scalar.mpfi", ["pydx/scalar/mpfi.c"], 
    include_dirs=[],
    libraries=["mpfi","mpfr","gmp"])

setup(
    name='pydx',
    version='0.4',
    description='Python Differential Geometry Tools',
    author='Simon Burton',
    author_email='simon@arrowtheory.com',
    url='http://gr.anu.edu.au/svn/people/sdburton/pydx',
    packages=['pydx', 'pydx.scalar', 'pydx.test'],
    ext_modules=[ext],
)


