#!/usr/bin/env python

#    manifold.py: (semi-)riemannian geometry.
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
from random import random

from pydx.scalar import Scalar
from pydx.mjet import MJet
from pydx.tensor import Tensor
from pydx.field import TensorField

__module__ = sys.modules[__name__]
Scalar.clients.append(__module__)

class RManifold(object):
    scalar_one = 1.0
    def __init__( self, g, g_uu=None ):
        if g_uu is not None:
            g.uu = g_uu
        g.uu.dd = g
        dim = g.dim
        g.p = g.comma()
        gamma = (self.scalar_one/2) * g.uu.mul( g.p + g.p.transpose(2,1,0) - g.p.transpose(2,0,1), (1,1) ).transpose(0,2,1)
        christoffel = gamma
        gamma.p = gamma.comma()
        # XX write a test for kretschman in SwartzchildMetric XX
        gamma_gamma = gamma.mul(gamma, (0,2))
        riemann = \
                gamma.p.transpose(0,2,3,1) \
            - gamma.p.transpose(0,2,1,3) \
            + gamma_gamma.transpose(2,0,3,1) \
            - gamma_gamma.transpose(2,0,1,3)
        ricci = riemann.contract( (0,2) )
        ricci.ud = ricci.mul( g.uu, (0,0) ) # ricci.ud
        curvature = ricci.ud.contract( (0,1) )
        riemann.dddd = riemann.mul( g, (0,0) ).transpose(3,0,1,2)
        riemann.uuuu = riemann.mul( g.uu, (1,0) ).transpose(0,3,1,2)
        riemann.uuuu = riemann.uuuu.mul( g.uu, (2,0) ).transpose(0,1,3,2)
        riemann.uuuu = riemann.uuuu.mul( g.uu, (3,0) )
        kretschman = riemann.uuuu.mul( riemann.dddd, (0,0), (1,1), (2,2), (3,3) )
        del gamma_gamma, g_uu
        self.__dict__.update( locals() )
Scalar.clients.append( RManifold )
