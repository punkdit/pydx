#!/usr/bin/env python

#    metric.py
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


from pydx.options import options
from pydx.mjet import MJet, Jet
from pydx.scalar import Scalar, set_symbolic_scalar, set_interval_scalar, restore_scalar
from pydx.scalar.symbolic import VarFloat
from pydx.scalar.mpfi import Interval, DisorderException
from pydx.tensor import Tensor
from pydx.field import TensorField
from pydx.scalar.fmath import exp, log, sin, cos, sqrt, sqr
from pydx.manifold import RManifold
from pydx.geodesic import Geodesic

#
###############################################################################
#

class HyperbolicMetric(TensorField):
    def __init__( self ):
        TensorField.__init__( self, (Tensor.dn,Tensor.dn), 2 )
        self.uu = self.Inverse()
        self.uu.dd = self
    def __call__( self, x, y ):
        zero = self.scalar_zero
        one = self.scalar_one
        g = [ [one/(y**2), zero], [zero, one/(y**2)] ]
        g = Tensor( self.valence, 2, elems=g )
        return g
    class Inverse(TensorField):
        def __init__( self ):
            TensorField.__init__( self, (Tensor.up,Tensor.up), 2 )
        def __call__( self, x, y ):
            zero = self.scalar_zero
            one = self.scalar_one
            g = [ [y**2, zero], [zero, y**2] ]
            g = Tensor( self.valence, 2, elems=g )
            return g
        
class SwartzchildMetric(TensorField):
    def __init__( self, M ):
        M = self.scalar_promote(M)
        self.M = M
        self.r = 2.0*M # outer event horizon
        self.uu = SwartzchildMetric.Inverse(M)
        TensorField.__init__( self, (Tensor.dn, Tensor.dn), 4 )
        self.uu.dd = self
    def __call__( self, t, r, theta, phi ):
#        print "metric", t, r, theta, phi
        zero = self.scalar_zero
#        print "metric zero:", zero
#        assert r > self.r, (r,self.r)
#        assert not theta.sin().contains(0.0), "not implemented..."
        M = self.M
        dt_dt = (2*M-r)/r
        dr_dr = r/(r-2*M)
        dtheta_dtheta = r**2
        dphi_dphi = (r*sin(theta))**2
        g = [
            [ dt_dt,   zero,   zero,          zero,      ],
            [ zero,    dr_dr,  zero,          zero,      ],
            [ zero,    zero,   dtheta_dtheta, zero,      ],
            [ zero,    zero,   zero,          dphi_dphi, ],
        ]
#        print 'g =', g
#        print
        g = Tensor( self.valence, self.dim, elems=g )
        return g
    class Inverse(TensorField):
        def __init__( self, M ):
            M = self.scalar_promote(M)
            self.M = M
            self.r = 2.0*M # outer event horizon
            TensorField.__init__( self, (Tensor.up, Tensor.up), 4 )
        def __call__( self, t, r, theta, phi ):
    #    print "metric", t, r, theta, phi
            zero = self.scalar_zero
            one = self.scalar_one
    #    print "metric zero:", zero
#            assert r > self.r, (r,self.r)
#            assert not theta.sin().contains(0.0), "not implemented..."
            M = self.M
            dt_dt = (2*M-r)/r
            dr_dr = r/(r-2*M)
            dtheta_dtheta = r**2
            dphi_dphi = (r*sin(theta))**2
            g = [
                [ one/dt_dt,   zero,   zero,          zero,      ],
                [ zero,    one/dr_dr,  zero,          zero,      ],
                [ zero,    zero,   one/dtheta_dtheta, zero,      ],
                [ zero,    zero,   zero,          one/dphi_dphi, ],
            ]
    #    print 'g =', g
    #    print
            g = Tensor( self.valence, self.dim, elems=g )
            return g

class KerrMetric(TensorField):
    def __init__( self, M, a ):
        M = self.scalar_promote(M)
        a = self.scalar_promote(a)
        self.M = M
        assert 0.0<abs(a)<=M # why not zero ??
        self.a = a
        self.r_plus = M + sqrt( sqr(M) - sqr(a) ) # outer event horizon
        self.uu = self.Inverse(M,a)
        TensorField.__init__( self, (Tensor.dn, Tensor.dn), 4 )
        self.uu.dd = self
    def __call__( self, t, r, theta, phi ):
#        assert r > self.r_plus, (r,self.r_plus)
#        assert not sin(theta).contains(0.0), "not implemented..."
        zero = self.scalar_zero
        a = self.a
        M = self.M
        Sigma = sqr(r) + sqr(a*cos(theta))
        Delta = sqr(r) + sqr(a) - 2.0*M*r
        A = sqr(sqr(r) + sqr(a)) - Delta * sqr( a * sin(theta) )
        dt_dt = -( 1.0 - 2.0*M*r / Sigma )
        dt_dphi = - ( 2.0 * M * r / Sigma ) * a * sqr(sin(theta))
        dr_dr = Sigma / Delta
        dtheta_dtheta = Sigma
        dphi_dphi = A * sqr(sin(theta)) / Sigma
        g = [
            [ dt_dt,   zero,   zero,          dt_dphi,   ],
            [ zero,    dr_dr,  zero,          zero,      ],
            [ zero,    zero,   dtheta_dtheta, zero,      ],
            [ dt_dphi, zero,   zero,          dphi_dphi, ],
        ]
#        print 'g =', g
#        print
        g = Tensor( self.valence, self.dim, elems=g )
        return g
    class Inverse(TensorField):
        def __init__( self, M, a ):
            M = self.scalar_promote(M)
            a = self.scalar_promote(a)
            self.M = M
            assert 0.0<abs(a)<=M # why not zero ??
            self.a = a
            self.r_plus = M + sqrt( sqr(M) - sqr(a) ) # outer event horizon
            TensorField.__init__( self, (Tensor.up, Tensor.up), 4 )
        def __call__( self, t, r, theta, phi ):
    #    assert r > self.r_plus, (r,self.r_plus)
    #    assert not sin(theta).contains(0.0), "not implemented..."
            zero = self.scalar_zero
            a = self.a
            M = self.M
            Sigma = sqr(r) + sqr(a*cos(theta))
            Delta = sqr(r) + sqr(a) - 2.0*M*r
            A = sqr(sqr(r) + sqr(a)) - Delta * sqr( a * sin(theta) )
            dt_dt = -( 1.0 - 2.0*M*r / Sigma )
            dt_dphi = - ( 2.0 * M * r / Sigma ) * a * sqr(sin(theta))
            dr_dr = Sigma / Delta
            dtheta_dtheta = Sigma
            dphi_dphi = A * sqr(sin(theta)) / Sigma
            g = [
                [ dt_dt,   zero,   zero,          dt_dphi,   ],
                [ zero,    dr_dr,  zero,          zero,      ],
                [ zero,    zero,   dtheta_dtheta, zero,      ],
                [ dt_dphi, zero,   zero,          dphi_dphi, ],
            ]
    #    print 'g =', g
    #    print
            Z = 1.0 / (g[0][0]*g[3][3] - g[0][3]*g[3][0])
            g_up = [
                [ g[3][3]*Z, zero, zero, -g[0][3]*Z, ],
                [ zero, 1.0/g[1][1], zero, zero, ],
                [ zero, zero, 1.0/g[2][2], zero, ],
                [ -g[3][0]*Z, zero, zero, g[0][0]*Z, ],
            ]
            g_up = Tensor( self.valence, self.dim, elems=g_up )
            return g_up

class SpatialCurzonMetric(TensorField):
    def __init__( self, m ):
        self.m = m
        self.uu = self.Inverse(m)
        TensorField.__init__( self, (Tensor.dn, Tensor.dn), 3 )
        self.uu.dd = self
    def __call__( self, r, z, phi ):
        zero = self.scalar_zero
        m = self.m
        R2 = sqr(r)+sqr(z)
        nu = -sqr(m)*sqr(r) / (2.0 * sqr(R2))
        lmda = -m/sqrt(R2)
        # dt_dt = -exp(2*lmda)
        c = exp(2*(nu-lmda))
        dr_dr = c
        dz_dz = c
        dphi_dphi = sqr(r)*exp(-2.0*lmda)
        g = [
            [ dr_dr,  zero,  zero,      ],
            [ zero,   dz_dz, zero,      ],
            [ zero,   zero,  dphi_dphi, ],
        ]
        g = Tensor( self.valence, self.dim, elems=g )
        return g
    class Inverse(TensorField):
        def __init__( self, m ):
            self.m = m
            TensorField.__init__( self, (Tensor.up, Tensor.up), 3 )
        def __call__( self, r, z, phi ):
            zero = self.scalar_zero
            m = self.m
            R2 = sqr(r)+sqr(z)
            nu = -sqr(m)*sqr(r) / (2.0 * sqr(R2))
            lmda = -m/sqrt(R2)
            # dt_dt = -exp(2*lmda)
            c = exp(2*(nu-lmda))
            dr_dr = c
            dz_dz = c
            dphi_dphi = sqr(r)*exp(-2.0*lmda)
            g = [
                [ dr_dr,  zero,  zero,      ],
                [ zero,   dz_dz, zero,      ],
                [ zero,   zero,  dphi_dphi, ],
            ]
    
            g_up = [
                [ 1.0/g[0][0], zero, zero, ],
                [ zero, 1.0/g[1][1], zero, ],
                [ zero, zero, 1.0/g[2][2], ],
            ]
            g_up = Tensor( self.valence, self.dim, elems=g_up )
            return g_up

class RZCurzonMetric(TensorField):
    # r and z coords only
    def __init__( self, m ):
        self.m = m
        self.uu = self.Inverse(m)
        TensorField.__init__( self, (Tensor.dn, Tensor.dn), 2 )
        self.uu.dd = self
    def __call__( self, r, z ):
#        print "RZCurzonMetric: r,z =", r, z
        zero = self.scalar_zero
        m = self.m
        R2 = sqr(r)+sqr(z)
        self.nu = -sqr(m)*sqr(r) / (2.0 * sqr(R2))
        r = R2[(0,)*R2.rank]
        if type(r) ==Interval:
            assert R2[(0,)*R2.rank]>0.0, R2[(0,)*R2.rank]
        self.lmda = -m/sqrt(R2)
        # dt_dt = -exp(2*lmda)
        c = exp(2*(self.nu-self.lmda))
        dr_dr = c
        dz_dz = c

        g = [
            [ dr_dr,  zero,  ],
            [ zero,   dz_dz, ],
        ]
#        print "RZCurzonMetric: c=", c.scalar_component
        g = Tensor( self.valence, self.dim, elems=g )
        return g
    class Inverse(TensorField):
        def __init__( self, m ):
            self.m = m
            TensorField.__init__( self, (Tensor.up, Tensor.up), 2 )
        def __call__( self, r, z ):
#            print "RZCurzonMetric.Inverse", r, z
            zero = self.scalar_zero
            one = self.scalar_one
            m = self.m
            R2 = sqr(r)+sqr(z)
            self.nu = -sqr(m)*sqr(r) / (2.0 * sqr(R2))
            self.lmda = -m/sqrt(R2)
            # dt_dt = -exp(2*lmda)
            c = exp(2*(self.nu-self.lmda))
            g_up = [
                [ one/c, zero, ],
                [ zero, one/c, ],
            ]
            g_up = Tensor( self.valence, self.dim, elems=g_up )
#            print "RZCurzonMetric.Inverse g_up[0,0]=", g_up[0,0]
            return g_up

