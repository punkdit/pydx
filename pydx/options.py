#!/usr/bin/env python

#    options.py
#    Copyright (C) 2006, Simon David Burton and The Australian Natonal University.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys
from random import seed
from scalar import Scalar

class Options(object):
    def __init__(self,**defaults):
        self.argv = sys.argv[:]
        for key,val in defaults.items():
            self.__dict__[key]=val
        for arg in sys.argv[1:]:
            if arg.count('='):
                lhs,rhs=arg.split('=')
                try:
                    self.__dict__[lhs]=eval(rhs)
                except NameError:
                    self.__dict__[lhs]=rhs # it's a string
            else:
                self.__dict__[arg]=True
    def __getattr__(self, name):
        return None
    def get(self, name, default=None):
        val = getattr(self,name)
        if val is None:
            val = default
        return val
    def dumpstate(self, file):
        for key,val in self.__dict__.items():
            if type(val)==__builtins__.file:
                continue
            file.write( '%s = %s\n' %( key, repr(val) ))

def main():
    global options

    seed(0)
    Scalar.set( float, float, 0.0, 1.0 )

    revision = 'Unknown'
    try:
        lines = os.popen('svn info').readlines()
        revision = [line for line in lines if line.count('Revision')][0]
        revision = revision.split(' ')[1].strip()
        revision = int(revision)
    except:
        pass

    options = Options(
        prec=64, step_size=0.0001, order=14, compile=True, steps=32,
        test='rz', verbose=3, resume=True, filename=None, revision=revision,
        _file=sys.stdout,
    )

    for name, testfunc in globals().items():
        #print name, testfunc, options.test
        if name.startswith('test_') and callable(testfunc) and name.count( options.test ):
            print "test=%s"%repr(name)
            options.testfunc=name
            if options.filename:
                if options.resume:
                    try:
                        file = open(options.filename,'r')
                    except IOError:
                        options.resume = False
                        file = open(options.filename,'w')
                else:
                    file = open(options.filename,'w')
            else:
                options.resume = False
                file = sys.stdout
            options._file = file
            if not options.resume:
                options.dumpstate(file)

            try:
                testfunc()
            except KeyboardInterrupt:
                pass
main()

#    test_jet_taylor_univariate()
#    test_jet_taylor_multivariate()
#    test_slice()
#    test_derive()
#    test_derive_multivariate()

#    test_tensor()
#    test_get_func()
#    test_rmanifold()

#    test_ode()
#    test_geodesic_snell()
#    test_geodesic_flatspace() # doesn't actually test anything

#    test_hyperbolic_geodesic()
#    test_swartzchild_geodesic()

#    test_rzcurzon_geodesic()

