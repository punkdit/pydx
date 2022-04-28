
from random import random

from pydx.scalar import set_mpf_scalar, set_interval_scalar, set_symbolic_scalar, restore_scalar
from pydx.scalar.mpfi import Interval
from pydx.mjet import MJet, Jet
from pydx.manifold import RManifold
from pydx.tensor import Tensor
from pydx.field import TensorField, ConcreteTensorField
from pydx.ode import ODE
from pydx.geodesic import Geodesic
from pydx.transform import MobiusTransform, LinearTransform

def test_geodesic_snell():
    # from "Reflections on Relativity": 
    # http://www.mathpages.com/rr/s8-04/8-04.htm

#    set_mpf_scalar()
    set_interval_scalar()

    scalar_promote = MJet.scalar_promote
    scalar_one = MJet.scalar_one
    scalrand = lambda : scalar_promote(random())

#    set_symbolic_scalar()
#    t_uu = Tensor( (Tensor.up,Tensor.up), 2 )
#    for idxs in t_uu.genidx():
#        t_uu[idxs] = VarFloat('t^%s' % (''.join(str(i) for i in idxs)))
#    t_ddd = Tensor( (Tensor.dn,Tensor.dn,Tensor.dn), 2 )
#    for idxs in t_ddd.genidx():
#        t_ddd[idxs] = VarFloat('t_%s' % (''.join(str(i) for i in idxs)))
#    gamma = t_uu.mul( t_ddd + t_ddd.transpose(2,1,0) - t_ddd.transpose(2,0,1), (1,1) ).transpose(0,2,1)
#    for idxs in gamma.genidx():
#        print 'gamma_%s' % (''.join(str(i) for i in idxs)), gamma[idxs]

#    return

#    restore_scalar()

    class ODE_for_func1(ODE):
        def __init__( self, x0, y0, A, B ):
            ODE.__init__( self, x0, y0 )
            self.A = A
            self.B = B
        def __call__( self, t, x ):
            assert len(x)==self.dim, x
            dx,dy = x[self.dim/2:] # these are the derivatives of x
            x,y = x[:self.dim/2]
            C = A/(A*x+B)
            result = Tensor( (Tensor.up,), self.dim )
            result[0] = dx
            result[1] = dy
            result[2] = -C * dx*dx + C * dy*dy
            result[3] = -2.0 * C * dx * dy
            return result

        
    A = scalar_promote(5.0)
    B = 1.0/A
    n_func0 = lambda x,y: scalar_one
    n_func1 = lambda x,y: A*x+B
    n_func2 = lambda x,y: scalar_one/y
    n_evals = 0
    for n_func in (n_func0, n_func1, n_func2):
            print '--'*30
            g = ConcreteTensorField.zero( (Tensor.dn,Tensor.dn), 2 )
            g.uu = ConcreteTensorField.zero( (Tensor.up,Tensor.up), 2 )
            for i in (0,1):
                g[i,i] = lambda x,y: n_func(x,y)*n_func(x,y)
                g.uu[i,i] = lambda x,y: scalar_one / (n_func(x,y)*n_func(x,y))
            manifold = RManifold( g )
            gamma = manifold.gamma
            
#            x, y = x0, y0 = (scalrand()-0.5)*10.0, (scalrand()-0.5)*10.0
            x, y = x0, y0 = scalar_promote(1.0), scalar_promote(1.0)
            dx, dy = dx0, dy0 = scalrand(), scalrand()-0.5
            print "x, y, dx, dy", x, y, dx, dy
            t = t0 = scalar_promote(0.0)

            if 1:
                n = n_func(x,y)
                gp_xy = g.p(x,y).apply( lambda r:r[()] )
                g00_func = lambda x,y: n_func(x,y)*n_func(x,y)
                gp_000 = MJet(1).promote( g00_func(Jet([x,scalar_one]), y) )[1]
                gp_001 = MJet(1).promote( g00_func(x, Jet([y,scalar_one])) )[1]
                assert abs( gp_xy[0,0,0] - gp_000 ) < 1e-10
                assert abs( gp_xy[0,0,1] - gp_001 ) < 1e-10
                assert abs( gp_xy[1,0,0] ) < 1e-10
                assert abs( gp_xy[0,1,0] ) < 1e-10
                assert abs( gp_xy[1,0,1] ) < 1e-10
                assert abs( gp_xy[0,1,1] ) < 1e-10
                assert abs( gp_xy[1,1,1] - gp_001 ) < 1e-10
                assert abs( gp_xy[1,1,0] - gp_000 ) < 1e-10
    
    #      for idxs in g.uu.genidx():
    #        print 'g^%s' % (''.join(str(i) for i in idxs)), g.uu(x,y)[idxs]
    #      for idxs in g.p.genidx():
    #        print 'g.p_%s' % (''.join(str(i) for i in idxs)), g.p(x,y)[idxs][()]
    #      for idxs in gamma.genidx():
    #        print 'gamma_%s' % (''.join(str(i) for i in idxs)), gamma(x,y)[idxs][()]
                gamma_000 = (scalar_one/n) * MJet(1).promote(n_func(Jet([x,scalar_one]),y))[1]
                gamma_001 = (scalar_one/n) * MJet(1).promote(n_func(x,Jet([y,scalar_one])))[1]
    
                assert abs( gamma(x,y)[0,0,0][()] - gamma_000 ) < 1e-10
                assert abs( gamma(x,y)[0,1,1][()] - -gamma_000 ) < 1e-10
                assert abs( gamma(x,y)[1,0,1][()] - gamma_000 ) < 1e-10
                assert abs( gamma(x,y)[1,1,0][()] - gamma_000 ) < 1e-10
                assert abs( gamma(x,y)[0,0,1][()] - gamma_001 ) < 1e-10
                assert abs( gamma(x,y)[0,1,0][()] - gamma_001 ) < 1e-10
                assert abs( gamma(x,y)[1,0,0][()] - -gamma_001 ) < 1e-10
                assert abs( gamma(x,y)[1,1,1][()] - gamma_001 ) < 1e-10

            STEPS = 10
            steps1 = []
            invarient = None

            ode = Geodesic( manifold, t0, [x,y,dx,dy] )
            h = scalar_promote(0.01)
            for i in range(STEPS):
#                print x, y, dx, dy
                n = n_func(x,y)
                if n_func in (n_func0,n_func1):
                    q = dy/(dx**2+dy**2).sqrt()
                else:
                    q = dx/(dx**2+dy**2).sqrt()
                _invarient = n*q
                print invarient, _invarient
                if invarient is None:
                    invarient = _invarient
                assert invarient.overlapping(_invarient)
                n_evals += 1

                x = Interval(x).lower
                Y0 = ode.contract1( Interval(x,(x+h).upper), ode.y )
                ode.istep2( h, Y0 )
                x, y, dx, dy = ode.y
                steps1.append( (x,y,dx,dy) )

            if n_func == n_func1:
                x, y = x0, y0
                dx, dy = dx0, dy0
                print "x, y, dx, dy", x, y, dx, dy
                t = t0 = scalar_promote(0.0)
                steps2 = []
    
                ode = ODE_for_func1( t0, [x,y,dx,dy], A, B )
                for i in range(STEPS):
                    x = Interval(x).lower
                    Y0 = ode.contract1( Interval(x,(x+h).upper), ode.y )
                    ode.istep2( h, Y0 )
                    x, y, dx, dy = ode.y
                    steps2.append( (x,y,dx,dy) )
                    n = n_func(x,y)
#                    print x, y, dx, dy
                    if n_func in (n_func0,n_func1):
                        q = dy/(dx**2+dy**2).sqrt()
                    else:
                        q = dx/(dx**2+dy**2).sqrt()
                    _invarient = n*q
                    print invarient, _invarient
                    assert invarient.overlapping(_invarient)
                    n_evals += 1
    
                for i in range(STEPS):
                    for x0, x1 in zip( steps1[i], steps2[i] ):
                        assert x0.overlapping(x1)
                        n_evals += 1

    restore_scalar()
    print "test_geodesic_snell: n_evals:", n_evals




def test_geodesic_flatspace():
#    set_mpf_scalar()
    set_interval_scalar()

    scalar_promote = MJet.scalar_promote
    scalar_one = MJet.scalar_one
    scalrand = lambda : scalar_promote(random())

    metric = TensorField.identity( (Tensor.dn,Tensor.dn), 2 )
    metric_up = TensorField.identity( (Tensor.up,Tensor.up), 2 )

    n_evals = 0
    for _ in range(1):
#        a, b, c, d = [(scalrand()-0.5,scalrand()-0.5) for _ in range(4)]
        a, b, c, d = ( (0.0,-1.0), (0.0,1.0), (1.0,0.0), (1.0,0.0), )
        for T in (
#            LinearTransform(*[scalrand()-0.5 for _ in range(4)]),
#            ComplexTransform(a, b),
            MobiusTransform(a, b, c, d),
        ):
            print T

            g = metric.transform(T)
            g_uu = metric_up.transform(T)

#            x, y = x0, y0 = (scalrand()-0.5)*10.0, (scalrand()-0.5)*10.0
            x, y = x0, y0 = T(0.0,0.0)
#            dx, dy = dx0, dy0 = scalrand(), scalrand()
            dx, dy = dx0, dy0 = (0.0,1.0)
            print x, y, dx, dy
            t = t0 = scalar_promote(0.0)

            manifold = RManifold( g, g_uu )
            ode = Geodesic( manifold, t0, [x,y,dx,dy] )
            h = scalar_promote(0.0001)
            for i in range(5):
                x = Interval(x).lower
#                print x+h
                Y0 = ode.contract1( Interval(x,(x+h).upper), ode.y )
                ode.istep1( h, Y0 )
                tensor = ConcreteTensorField( (Tensor.dn,), 2 )
                tensor[0] = lambda x,y: dx
                tensor[1] = lambda x,y: dy
                x, y, dx, dy = ode.y

                print x, y

                for _tensor in (
#                    tensor(x,y),
#                    tensor.transform( T )(x,y),
#                    tensor.transform( T )(*T(x,y)),
#                    tensor.transform( T )(*T.inverse(x,y)),
#                    tensor.transform( T.inverse )(x,y),
#                    tensor.transform( T.inverse )(*T(x,y)),
#                    tensor.transform( T.inverse )(*T.inverse(x,y)),
                 ):
                    print MJet(0).promote(_tensor[0])[()],
                print

                n_evals += 1

#                print '\t', (dx, dy)
#                print "step"

    restore_scalar()
    print "test_geodesic_flatspace: n_evals:", n_evals


        
