
from random import random

from gmpy import mpf

from pydx.mjet import MJet, Jet, Derivation, multi_unit, multi_range_order, cross
from pydx.scalar import set_interval_scalar, set_symbolic_scalar, restore_scalar, Scalar
from pydx.scalar.mpfi import Interval
from pydx.scalar.symbolic import VarFloat

def test_jet_taylor_univariate():
    fs = []
    fs += [
        lambda x: x,
        lambda x: 1.0/(x+1.0),
        lambda x: x*x-1.0,
        lambda x: x**2-1.0,
#        lambda x: (x+1.0)**0.5, # does not work on float scalars
        lambda x: (x+1.0)*(x-2.0)*x,
        lambda x: (x+1.0)*(x-2.0)/(x+10.0),
        lambda x: (x+1.0)*(x+2.0)*x*(x+1.0)*(x+2.0)*x,
        lambda x: 7.0*x**3-2.0*x**2+1.0,
    ]

    n_evals = 0
         
    total_err = 0.0
    for f in fs:
        for trials in range(50):
            x = random()-0.5
            h = 0.01*random()
            order = 5
            dx = Jet([x,1.0])
            y = f(dx)
            r0 = f(x+h)
            r1 = y.expand( (h,), order )
            if abs(r0)+abs(r1)>1e-6:
                n_evals += 1
                err = abs(r0-r1)/(abs(r0)+abs(r1))
                assert err < 1e-6, (r0,r1)
                total_err += err
    print "total_err:", total_err

# These tests need interval scalar's to work:
    fs.append( lambda x: x.exp() )
    fs.append( lambda x: (2.0*x).exp() )
    fs.append( lambda x: (x**3).exp() )
    fs.append( lambda x: (x+1.0).log() )
    fs.append( lambda x: (x**3+1.0).log() )
    fs.append( lambda x: (x+1.0)**0.5 )
    fs.append( lambda x: (x**3+1.0).sin() )
    fs.append( lambda x: (x**3+1.0).cos() )
    fs.append( lambda x: (x**2).tan() )
    fs.append( lambda x: (x**2).asin() )
    fs.append( lambda x: (x**3+1.0).atan() )
    fs.append( lambda x: (x**3+1.0).sqrt() )
    fs.append( lambda x: (x**3+1.0).sinh() )
    fs.append( lambda x: (x**3+1.0).cosh() )
    fs.append( lambda x: (x**3+1.0).tanh() )
         
    set_interval_scalar()
    #assert len(Scalar.stack)==1

    for trials in range(4):
        for f in fs:
            x = Interval(random()-0.5)
            h = Interval(random()-0.5)
            xx = x.hull(x+h)
            for order in range(4):
                dx = Jet([x,1.0])
                dxx = Jet([xx,1.0])
                err = f(dxx)[order+1] 
                y = f(dx)
                r0 = f(x+h)
                r1 = y.expand_err( (h,), order, (err,) )
                assert r1.overlapping(r0)
                n_evals += 1

    restore_scalar()

    print "test_jet_taylor_univariate n_evals:", n_evals
        
def test_jet_taylor_multivariate():
    fs = [
#        lambda x: x,
#        lambda x: x*x+1.0,
#        lambda x: x**2-1.0,
        lambda x: (x+1.0)*(x+-2.0)*x,
#        lambda x: (x+1.0)*(x-2.0)/(x*x+10.0),
#        lambda x: (x+1.0)*(x+-2.0)*x*(x+0.5),
#        lambda x: 5.0*x**3-2.0*x**2+1.0,
    ]
    gs =  []
#    gs += [ lambda x, y, f1=f1, f2=f2: f1(x) for f1 in fs for f2 in fs ]
    gs += [ lambda x, y, f1=f1, f2=f2: f1(x)+f2(y) for f1 in fs for f2 in fs ]
    gs += [ lambda x, y, f1=f1, f2=f2: f1(x)*f2(y) for f1 in fs for f2 in fs ]
#    gs += [ lambda x, y, f1=f1, f2=f2: f1(x)+f2(y)*f1(y) for f1 in fs for f2 in fs ]
#    gs += [ lambda x, y, f1=f1, f2=f2: f1(x)*f2(y)+f1(y) for f1 in fs for f2 in fs ]

    max_err = 0.0
    n_evals = 0

    for g in gs:
        for trials in range(30):
            x = random()-0.5, random()-0.5
            h = (0.01*random(),0.01*random())
            dx = Jet({ (0,0):x[0], (1,0):1.0 }), Jet({ (0,0):x[1], (0,1):1.0 })
            y = g(*dx)
            r0 = g( x[0]+h[0], x[1]+h[1] )
#            r2 = y.expand( h, 2 )
            r3 = y.expand( h, 3 )
            r = abs(r0)+abs(r3)
            if r>1e-6:
                err = abs(r0-r3)/r 
                assert err < 1e-5, err
                n_evals += 1
    print "test_jet_taylor_multivariate n_evals:",n_evals

    # The following ops are only defined on intervals
    for g in gs[:]:
        gs.append( lambda x, y, g=g: g(x,y).exp() )
        gs.append( lambda x, y, g=g: (g(x,y)**2+0.1).log() )
        gs.append( lambda x, y, g=g: (g(x,y)**2+0.1)**1.1 )
        gs.append( lambda x, y, g=g: g(x,y).sin() )
        gs.append( lambda x, y, g=g: g(x,y).cos() )
        gs.append( lambda x, y, g=g: g(x,y).tan() )
##     gs.append( lambda x, y, g=g: g(x,y).acos() ) # not implemented
        gs.append( lambda x, y, g=g: g(x,y).atan() )
        gs.append( lambda x, y, g=g: g(x,y).sinh() )
        gs.append( lambda x, y, g=g: g(x,y).cosh() )
        gs.append( lambda x, y, g=g: g(x,y).tanh() )

        gs.append( lambda x, y, g=g: g(x,y)**3 )

    gs.append( lambda x, y, g=g: ((x+y)/2.0).asin() )
    gs.append( lambda x, y, g=g: (2.0+(x+y)).sqrt() )

    set_interval_scalar(128)

    total_width = mpf(0.0)
    for trials in range(30):
        for g in gs:
            x = [ Interval(random()-0.5) for i in (0,1) ]
            h = [ Interval(random()-0.5)*0.01 for i in (0,1) ]
            xx = [ x[i].hull(x[i]+h[i]) for i in (0,1) ]
            for order in (0,1,2):
                r0 = g(x[0]+h[0],x[1]+h[1])
                dx = [ Jet( {(0,0):x[i],multi_unit(2,i):1.0} ) for i in (0,1) ]
                y = g(*dx)
                dxx = [ Jet( {(0,0):xx[i],multi_unit(2,i):1.0} ) for i in (0,1) ]
                yy = g(*dxx)
                err = [ yy[j] for j in multi_range_order(2,order+1,order+1) ]
                r1 = y.expand( h, order )
                r1 = y.expand_err( h, order, err )
                total_width += r1.width()
                assert r1.overlapping(r0), (x,h)
                n_evals += 1
    print "total_width:", total_width

    restore_scalar()

    print "test_jet_taylor_multivariate n_evals:",n_evals
        

#
###############################################################################
#

def test_slice():
    r = 1.234
    x = Jet( {(0,0):r, (1,0):1.0} )
    x0 = Jet( {(0,):r, (1,):1.0} )
    x1 = Jet( {(0,):0.0, (1,):0.0} )
    assert x[:,0] == x0
    assert x[:,1] == x1
    assert x0 != x1

def derive( f, *y ):
    return Derivation(f)(*y)

def test_derive():

    tests = [
    # list of (f, f')
        ( lambda x: x, lambda x: 1.0 ),
        ( lambda x: x*x, lambda x: 2.0*x ),
        ( lambda x: x*x+7.0, lambda x: 2.0*x ),
        ( lambda x: 3.0*x*x-9.0*x, lambda x: 6.0*x-9.0 ),
        ( lambda x: 2.0*x**3+1.0, lambda x: 6.0*x*x ),
    ]

    n_evals = 0
    for f, df in tests:
        for i in range(5):
            # scalar
            x = random()-0.5
            y0 = derive( f, x )[0]
            y0 = y0[()] # XX should this be a scalar already ?
            y1 = df( x )
            assert abs(y0-y1)<1e-10, abs(y0-y1)
            n_evals += 1
    
            # rank 1
            rank = 1
            idxs = (0,)
            x = Jet( {(0,):x,(1,):1.0} )
            y0 = derive( f, x )[0]
            y0 = y0[idxs]
            y1 = df( x )
            y1 = MJet(rank).promote(y1)
            y1 = y1[idxs]
            assert abs(y0-y1)<1e-10, abs(y0-y1)
            n_evals += 1
 
            # rank 1 to 4
            for rank in range(1,5):
                x = Jet( rank = rank )
                for idxs in cross( (3,)*rank ):
                    x[idxs] = random()-0.5
                y0 = derive( f, x )[0]
                y1 = df( x )
                y1 = MJet(rank).promote(y1)
                for idxs in cross( (3,)*rank ):
                    _y0 = y0[idxs]
                    _y1 = y1[idxs]
                    assert abs( _y0 - _y1 ) < 1e-10, abs(_y0-_y1)
                    n_evals += 1
    print "test_derive: n_evals:", n_evals

def test_derive_multivariate():
    _tests = [
    # list of (f, f')
        ( lambda x: x, lambda x: 1.0 ),
#        ( lambda x: x*x, lambda x: 2.0*x ),
        ( lambda x: x*x+7.0, lambda x: 2.0*x ),
        ( lambda x: 3.0*x*x-9.0*x, lambda x: 6.0*x-9.0 ),
        ( lambda x: 2.0*x**3+1.0, lambda x: 6.0*x*x ),
    ]
    tests = []
    for f1, df1 in _tests:
        for f2, df2 in _tests:
#            tests.append(
#                ( lambda x, y: f1(x)+f2(x), lambda x, y: df1(x)+df2(x), lambda x, y: 0.0 ))
#            tests.append(
#                ( lambda x, y: f1(y)+f2(y), lambda x, y: 0.0, lambda x, y: df1(y)+df2(y) ))
            tests.append(
                ( lambda x, y: f1(x)+f2(y), lambda x, y: df1(x), lambda x, y: df2(y) ))
            tests.append(
                ( lambda x, y: (f1(x)+f2(y))**2, 
                    lambda x, y: 2.0*(f1(x)+f2(y))*df1(x), 
                    lambda x, y: 2.0*(f1(x)+f2(y))*df2(y), ))
#            tests.append(
#                ( lambda x, y: f1(x)*f2(y), lambda x, y: df1(x)*f2(y), lambda x, y: f1(x)*df2(y) ))
    n_evals = 0
    for f, fx, fy in tests:
#        for _ in range(3):
            # test scalar args
            x = [ random()-0.5 for _ in (0,1) ]
            df = derive( f, *x )
            for i in (0,1):
                y0 = df[i][()]
                y1 = (fx,fy)[i]( *x )
                assert abs(y0-y1)<1e-10
                n_evals += 1
    
            # test rank-1 arg
            x = [ random()-0.5 for _ in (0,1) ]
            rank = 1
            idxs = (0,)
            x = [ Jet( {(0,):xi,(1,):1.0} ) for xi in x ]
            df = derive( f, *x )

            for i in (0,1):
                y0 = df[i][(0,)]
                y1 = (fx,fy)[i]( *x )
                y1 = MJet(rank).promote(y1)
                y1 = y1[idxs]
                assert abs(y0-y1)<1e-10
                n_evals += 1
    
            # test rank-n args
            for rank in range(1,3):
                x = [ Jet( rank = rank ) for _ in (0,1) ]
                for xi in x:
                    for idxs in cross( (3,)*rank ):
                        xi[idxs] = random()-0.5
                y0 = derive( f, *x )
                for i in (0,1):
                    y1 = (fx,fy)[i]( *x )
                    y1 = MJet(rank).promote(y1)
                    for idxs in cross( (3,)*rank ):
                        _y0 = y0[i][idxs]
                        _y1 = y1[idxs]
                        assert abs( _y0 - _y1 ) < 1e-10
                        n_evals += 1
 
    print "test_derive_multivariate: n_evals:", n_evals
