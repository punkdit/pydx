#!/usr/bin/env python

#    ode.py
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


from pydx.mjet import MJet, Jet, Spray
from pydx.scalar import Scalar, set_symbolic_scalar, set_interval_scalar, restore_scalar
from pydx.scalar.symbolic import VarFloat, ComputationContext
from pydx.scalar.mpfi import Interval, DisorderException
from pydx.tensor import Tensor
from pydx.field import TensorField

#
###############################################################################
#


class Info(object):
    def __init__(self,**kw):
        self.__dict__.update(globals())
        self.__dict__.update(kw)
    def __getitem__(self, name):
        return self.__dict__[name]

class ODE(object):
    def __init__( self, x0, y0, paramname='t', varnames=None ):
        # XX use get/setattr trickery to hide/expose the param/var attrs XX#
        self.x = x0
        self.y = y0
        ####################################################################
        self.dim = len(self.y)
        if varnames is None:
            varnames = ['x%d'%i for i in range(self.dim)]
        self.varnames = varnames
        self.paramname = paramname
        self.i = 0 # Step number
#        self.spray_funcs = {} # map order -> func
        self.step_funcs = {} # map order -> func
    def dumpstate(self, file):
        file.write( 'i = %s\n' % self.i )
        file.write( '%s = %s\n' % ( self.paramname, repr(self.x) ) )
        for idx, name in enumerate(self.varnames):
            file.write( '%s = %s\n' % ( name, repr(self.y[idx]) ) )
        file.flush()
    def loadstate(self, file):
        info = Info(i=self.i)
        while self.i==info.i:
            # read lines until the counter increments
            line = file.readline().strip()
            if not line:
                break
            exec line in info.__dict__
        try:
            y = self.y[:]
            x = info[self.paramname]
            for idx, name in enumerate(self.varnames):
                y[idx] = info[name]
            i = info.i
        except Exception, e:
            raise
        self.y[:] = y
        self.x = x
        self.i = i
    def __call__( self, x, y ):
        """
            ODE(x,y)
                @param x : parameter (scalar)
                @param y : vector

            y'(x) = f(x,y)
        """
        pass
    def step1( self, h ):
        h = self.scalar_promote(h)
        v = self.y
        t = self.x
        dv = self( t, v )
        for i, dvi in enumerate(dv):
            if isinstance(dvi,MJet):
                dv[i] = dvi[()]
        v = [ v[i] + h*dv[i] for i in range(self.dim) ]
        t = t + h
        self.y = v
        self.x = t
        self.i += 1
    def step2( self, h ):
        h = self.scalar_promote(h)
        v = self.y
        t = self.x
        dv = self( t, v )
        for i, dvi in enumerate(dv):
            if isinstance(dvi,MJet):
                dv[i] = dvi[()]
        v = [ Jet([v[i],dv[i]]) for i in range(self.dim) ]
        ddv = self( Jet( [t, self.scalar_one] ), v )
        ddv = [ (self.scalar_one/2) * ddvi[1] for ddvi in ddv ]
        hh = h*h
        v = [ v[i] + h*dv[i] + hh*ddv[i] for i in range(self.dim) ]
        t = t + h
        self.y = v
        self.x = t
        self.i += 1
    def istep1( self, h, yy ):
        """
            First order interval step.
            @param yy: bound on solution (got from contract1)
        """
        f = self
        h = self.scalar_promote(h)
        v = self.y
        t = self.scalar_promote( self.x )
        dv = f( t, yy )
        for i, dvi in enumerate(dv):
            if isinstance(dvi,MJet):
                dv[i] = dvi[()]
        v = [ v[i] + h*dv[i] for i in range(self.dim) ]
#            print "istep1: err:", (h*dv[0]).width()
        t = t + h
        self.y = v
        self.x = t
        self.i += 1
    def istep2( self, h, y ):
        """
            Second order interval step.
            @param y: bound on solution (got from contract1)
        """
        f = self
        h = self.scalar_promote(h)

        x0 = self.x
        y0 = self.y

        dy0 = f( x0, y0 )
        dy0 = [ MJet(0).promote(dy0i)[()] for dy0i in dy0 ]

        x = x0.hull(x0+h)
        dy = f( x, y )
        dy = [ MJet(0).promote(dyi)[()] for dyi in dy ]

        y_jet = [ Jet([y[i],dy[i]]) for i in range(self.dim) ]
        x_jet = Jet( [x, self.scalar_one] )

        ddy = f( x_jet, y_jet )
        ddy = [ (self.scalar_one/2) * ddyi[1] for ddyi in ddy ]

        hh = h**2
        y = [ y0[i] + h*dy0[i] + hh*ddy[i] for i in range(self.dim) ]
        x1 = x0 + h.lower
#        print "istep2: err:", (hh*ddy[0]).width()

        self.y = y
        self.x = x1
        self.i += 1

    def _compile_spray(self, order):
        f = self
        dumpfile = open('stepdump.py','w')
        set_symbolic_scalar()
        x = VarFloat('x')
        y = [ VarFloat( 'y%d'%i ) for i in range(self.dim) ]
        vs = [ x ] + y
        name = 'spray'
        args = vs + [ 'rval' ]
        ctx = ComputationContext(name, args, dumpfile)
        rval = [None]*self.dim
        lazy = Spray(f, x, y)
        for i in range(order):
            ctx.assign( 'rval[%d]'%i, lazy[i] )
        func = ctx.finalize()
        restore_scalar()
        return func
    def get_spray(self, order):
        if order in self.spray_funcs:
            return self.spray_funcs[order]
        func = self._compile_spray(order)
        self.spray_funcs[order]=func
        return self.spray_funcs[order]
    def compile(self,n):
        assert n not in self.step_funcs
        dumpfile = open('stepdump.py','w')
        set_symbolic_scalar()
        x0 = VarFloat('x0')
        y0 = [ VarFloat( 'y0_%d'%i ) for i in range(self.dim) ]
        h = VarFloat('h')
        x = VarFloat('x')
        y = [ VarFloat( 'y_%d'%i ) for i in range(self.dim) ]
        vs = [x0]+y0+[h,x]+y
        name = 'spray'
        args = vs + [ 'y1' ]
        print "compile n=",n,"args:",args
        f = self
        ################## START COMPUTATION
        print "computing..."
        ctx = ComputationContext(name, args, dumpfile)
        y0 = Spray( f, x0, y0 ) # a list of MJet's
        y = Spray( f, x, y ) # a list of MJet's
        y1 = [ y0i.expand_err( (h,), n-1, (yi[n],) ) for y0i,yi in zip(y0,y) ]
        for i in range(self.dim):
            ctx.assign( 'y1[%d]'%i, y1[i] )
        print "finalize..."
        func = ctx.finalize()
        ################## FINISH COMPUTATION
        restore_scalar()
        return func

    def istepn( self, h, y, n=2 ):
        """
            Arbitrary order interval step.
            @param h: step size
            @param y: bound on solution (got from contract1)
        """
#        print type(y)
#        print "istepn h=%s y=%s n=%d" % ( h, y, n )
        func = self.step_funcs.get(n)
        if func is None:
            func = self.compile(n)
            self.step_funcs[n] = func
        x0 = self.x
        y0 = self.y
        h = self.scalar_promote(h)
        x = x0.hull(x0+h)
        y1=[None]*self.dim
        args = [x0]+y0+[h,x]+[yi for yi in y]+[y1]
        func( *args )
        x1 = x0 + h

        self.y = y1
        self.x = x1
        self.i += 1
    def _istepn( self, h, y, n=2 ):
        """
            Arbitrary order interval step.
            @param h: step size
            @param y: bound on solution (got from contract1)
        """
        f = self
        h = self.scalar_promote(h)

        x0 = self.x
        y0 = self.y
        x = x0.hull(x0+h)

        y0 = Spray( f, x0, y0 ) # a list of MJet's
        y = Spray( f, x, y ) # a list of MJet's

        x1 = x0 + h
        y1 = [ y0i.expand_err( (h,), n-1, (yi[n],) ) for y0i,yi in zip(y0,y) ]

        self.y = y1
        self.x = x1
        self.i += 1
    #####################
    def contract1( self, x, y0, max_iter=100 ):
        """
        """
        set_interval_scalar()
        f = self
        hh = Interval(0.0, x.width())
        _y0 = Tensor( (Tensor.up,), len(y0) )
#        print "contract1", x, y0
        for i in range(len(_y0)):
            _y0[i] = ODE.scalar_promote(y0[i])
        _y = y = y0 = _y0

        cookies = 1 # makes no difference...
        
        i = 0
        while 1:
            y = _y
#            print "hh =",hh
#            for i in range(self.dim):
#                print "y0[%d] = %s"%(i,MJet(0).promote(y0[i])[()])
#                print "_y[%d] = %s"%(i,MJet(0).promote(_y[i])[()])
            try:
#                print '--'*30
#                print "y", y
#                print "f(x,y)",f(x,y)
#                print "hh*f(x,y)",hh*f(x,y)
                _y = y0 + hh * f(x,y) # Contract
            except DisorderException:
                raise
    #    print "y =", y
    #    print "_y =", _y
            shrink = True
            for i in range(self.dim):
                _y[i] = MJet(0).promote(_y[i])[()]
                if not y[i].contains(_y[i]):
                    shrink = False
            if shrink:
                cookies -= 1
            if cookies <= 0:
                break
            i += 1
            if i > max_iter:
                return None
#        print "contract1: contraction:", y
        restore_scalar()
        return y
Scalar.add_client( ODE )

