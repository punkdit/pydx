#!/usr/bin/env python

#    geodesic.py: the geodesic equation.
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

from pydx.mjet import MJet, Jet
from pydx.scalar import Scalar, set_symbolic_scalar, set_interval_scalar, restore_scalar
from pydx.scalar.mpfi import Interval

from pydx.tensor import Tensor
from pydx.field import TensorField, ConcreteTensorField
from pydx.manifold import RManifold
from pydx.transform import MobiusTransform
from pydx.options import options
from pydx.ode import ODE


class Geodesic(ODE):
    def __init__( self, manifold, x0, y0, paramname=None, varnames=None ):
        """
        """
        ODE.__init__( self, x0, y0, paramname=paramname, varnames=varnames )
        self.manifold = manifold
        self.gamma = manifold.gamma
        if options.compile:
            self.gamma = self.gamma.compile()
    def __call__( self, t, x ):
        """
            build Vector of derivatives of each component of Vector x
        """
        assert x is not None
        assert len(x)==self.dim, x
        dx = x[self.dim/2:] # these are the derivatives of x
        x = x[:self.dim/2]
        assert len(x)==self.dim/2
        assert len(dx)==self.dim/2
        result = Tensor( (Tensor.up,), self.dim )
        gamma = self.gamma( *x )
        for i in range(self.dim/2):
            result[i] = dx[i]
        for i in range(self.dim/2):
            r = sum( [ -gamma[i,j,k]*dx[j]*dx[k]
                        for j in range(self.dim/2) 
                        for k in range(self.dim/2) ], self.scalar_zero )
            result[self.dim/2+i] = r
        return result

