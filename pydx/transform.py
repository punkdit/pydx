#!/usr/bin/env python

#    transform.py
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


from random import random
from time import time

from pydx import scalar
from pydx.scalar import set_symbolic_scalar, restore_scalar
from pydx.scalar.symbolic import VarFloat
from pydx.tensor import Tensor, up, dn
from pydx.field import TensorField, ScalarField, Transform


#
###############################################################################
#
# Two dimensional Transforms follow..

class Complex(object):
    scalar_zero = 0.0
    scalar_promote = float
    def __init__( self, a, b=scalar_zero ):
        self.a = a # could be MJet, etc.
        self.b = b
    def __getitem__( self, idx ):
        return (self.a, self.b)[idx]
    def promote( cls, x ):
        if isinstance(x,Complex):
            return x
        elif isinstance(x,tuple):
            return Complex(*x)
        return Complex(x)
    promote = classmethod(promote)
    def __add__( self, other ):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Complex( a+c, b+d )
    def __sub__( self, other ):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Complex( a-c, b-d )
    def __neg__( self ):
        a, b = self.a, self.b
        return Complex( -a, -b )
    def __mul__( self, other ):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Complex( a*c-b*d, b*c+a*d )
    def reciprocal( self ):
        a, b = self.a, self.b
        norm = a*a+b*b
        return Complex( a/norm, -b/norm )
    def __div__( self, other ):
        return self * other.reciprocal()
scalar.Scalar.clients.append(Complex)

class LinearTransform(Transform):
    def __init__( self, a00, a01, a10, a11, inverse=None ):
        self.dim = 2
        self.a00 = a00
        self.a01 = a01
        self.a10 = a10
        self.a11 = a11
        if inverse is None:
            det = a00*a11 - a01*a10
            a00, a01, a10, a11 = a11/det, -a01/det, -a10/det, a00/det
            inverse = LinearTransform( a00, a01, a10, a11, self )
        assert inverse
        Transform.__init__( self, [self.call0,self.call1], inverse )
    def call0( self, x0, x1 ):
        return self.a00*x0 + self.a01*x1
    def call1( self, x0, x1 ):
        return self.a10*x0 + self.a11*x1

class ComplexTransform(Transform):
    def __init__( self, a, b, inverse=None ):
        self.dim = 2
        a = Complex.promote(a)
        b = Complex.promote(b)
        self.a = a
        self.b = b
        if inverse is None:
            a, b = Complex(1.0)/a, -b/a
            inverse = ComplexTransform( a, b, self )
        assert inverse
        Transform.__init__( self, [self.call0,self.call1], inverse )
    def call0( self, z0, z1 ):
        z = Complex(z0,z1)
        return (self.a*z+self.b)[0]
    def call1( self, z0, z1 ):
        z = Complex(z0,z1)
        return (self.a*z+self.b)[1]

class MobiusTransform(Transform):
    def __init__( self, a, b, c, d, inverse=None ):
        self.dim = 2
        self.a = Complex.promote(a)
        self.b = Complex.promote(b)
        self.c = Complex.promote(c)
        self.d = Complex.promote(d)
        assert self.a*self.d-self.b*self.c != Complex.promote(self.scalar_zero)
        if inverse is None:
            a, b, c, d = -self.d, self.b, self.c, -self.a
            inverse = MobiusTransform( a, b, c, d, self )
        Transform.__init__( self, [self.call0,self.call1], inverse )
    def call0( self, z0, z1 ):
        z = Complex(z0,z1)
        w = (self.a*z+self.b)/(self.c*z+self.d)
        return w[0]
    def call1( self, z0, z1 ):
        z = Complex(z0,z1)
        w = (self.a*z+self.b)/(self.c*z+self.d)
        return w[1]

