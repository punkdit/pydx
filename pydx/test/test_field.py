
from random import random
from time import time

from pydx import scalar
from pydx.scalar import set_symbolic_scalar, restore_scalar
from pydx.scalar.symbolic import VarFloat
from pydx.tensor import Tensor, up, dn
from pydx.field import TensorField, ScalarField, Transform
from pydx.transform import Complex, LinearTransform, ComplexTransform, MobiusTransform

#
###############################################################################
#

def test_tensor():
    t = Tensor( (dn,), 4 )
    t[0,] = 1.0
#    print t

    n_evals = 0
    for _ in range(10):
        a, b, c, d = [(random()-0.5,random()-0.5) for _ in range(4)]
        for T in (
            LinearTransform(*[random()-0.5 for _ in range(4)]),
            ComplexTransform(a, b),
            MobiusTransform(a, b, c, d),
        ):
        
            # test the inverse of the transform
            x, y = (random()-0.5)*10.0, (random()-0.5)*10.0
            _x, _y = T( x, y )
            _x, _y = T.inverse( _x, _y )
            assert abs(x-_x)+abs(y-_y) < 1e-6
            n_evals += 1

            # test .view
            _x, _y = T.inverse.view( T )( x, y )
            assert abs(x-_x)+abs(y-_y) < 1e-6
            n_evals += 1
            
            # test the inverse of the partial matrix
            # is the partial of the inverse transform
            Tx, Ty = T(x, y)
            p1 = T.partial(x,y)  # 2-2 matrix of partials
            p2 = T.inverse.partial(Tx, Ty)  # 2-2 matrix of partials of the inverse transform
            p = p1.mul( p2, (1,0) ) # sum over second index of p1 and first index of p2
            p.apply( lambda x:x[()] ) # get the scalar values out of the 0-Jet's
            # test identity matrix
            for r in p[1,1], p[1,1]:
                assert abs(r-1)<1e-10
            for r in p[0,1], p[1,0]:
                assert abs(r)<1e-10
            n_evals += 1

            # test ScalarField's are invarient
            a, b, c = [random()-0.5 for _ in range(3)]
            f = lambda x,y: a*x*x+b*y+c
            sf = ScalarField( 2, f )
            s1 = sf( *T(x,y) )
            s2 = sf.transform(T)( x, y )
            assert abs((s1-s2)[()])<1e-6
            n_evals += 1

            # test vector fields transform back 
            vf0 = sf.comma()
            vf1 = vf0.transform(T)
            vf2 = vf1.transform(T.inverse)
            assert vf0( x, y ).is_close( vf2(x, y) )
            n_evals += 1

            # test some tensor fields transform back 
            for valence in [ 
                    (Tensor.up,Tensor.up), 
                    (Tensor.up,Tensor.dn), 
                    (Tensor.dn,Tensor.dn), 
                    (Tensor.dn,Tensor.up), ]:
                k0 = TensorField.identity( valence, 2 )
                k1 = k0.transform(T)
                k2 = k1.transform(T.inverse)
                assert k0.valence == k1.valence == k2.valence
                assert k0( x, y ).is_close( k2(x, y) )
                n_evals += 1

            # test that vector fields transform correctly
            sf0 = sf
            vf0 = sf0.comma()
            sf1 = sf0.transform(T)
            vf10 = sf1.comma()
            vf01 = vf0.transform(T)
            assert vf10( x, y ).is_close( vf01(x, y) )
            n_evals += 1



#        print v(1.2), type(v(1.2)), v(1.2).valence
    print "test_tensor: n_evals:", n_evals


def test_get_func():
    scalar_zero = 0.0
    scalar_one = 1.0
    scalar_promote = float

    xys = [(random()-0.5,random()-0.5) for _ in range(100)]

    metric = TensorField.identity( (Tensor.dn,Tensor.dn), 2 )
    metric_up = TensorField.identity( (Tensor.up,Tensor.up), 2 )

    n_evals = 0
    T=MobiusTransform(*[(scalar_promote(1.2+i),scalar_promote(3.3*i)) for i in (-2,3,1,7)])
    g = metric.transform(T)
    g_uu = metric_up.transform(T)
    gp = g.comma()
    gamma = (1.0/2)*g_uu.mul( gp + gp.transpose(0,2,1) - gp.transpose(1,2,0), (1,0) )

    rs = []
    t0 = time()
    for x, y in xys:
#        r = g(x,y)[0,0][()]
        r = gamma(x,y)[0,0,0][()]
        rs.append(r)
    print "time:", time()-t0

    set_symbolic_scalar()

    metric = TensorField.identity( (Tensor.dn,Tensor.dn), 2 )
    metric_up = TensorField.identity( (Tensor.up,Tensor.up), 2 )

    n_evals = 0
    T=MobiusTransform(*[(scalar_promote(1.2+i),scalar_promote(3.3*i)) for i in (-2,3,1,7)])
    g = metric.transform(T)
    g_uu = metric_up.transform(T)
    gp = g.comma()
    gamma = (1.0/2)*g_uu.mul( gp + gp.transpose(0,2,1) - gp.transpose(1,2,0), (1,0) )
#    expr = g(VarFloat('x'),VarFloat('y'))[0,0][()]
    expr = gamma(VarFloat('x'),VarFloat('y'))[0,0,0][()]
    print "deeplen:", expr.deeplen()
    print "uniqlen:", expr.uniqlen()
    s_expr = expr.expr()
    print "strlen:", len(s_expr)
#    print s_expr

#    this breaks down on big expressions
#    t0 = time()
#    for i, (x, y) in enumerate(xys):
#        x = scalar_promote(x)
#        y = scalar_promote(y)
#        r = eval(s_expr, {'x':eval(str(x)),'y':eval(str(y))})
#        assert abs(r-rs[i])<1e-10
#    print "time:", time()-t0

    func = expr.get_func( 'func', ['x','y'] )
    t0 = time()
    for i, (x, y) in enumerate(xys):
        r = func(x,y)
        assert abs(r-rs[i])<1e-10
    print "time:", time()-t0

    restore_scalar()

