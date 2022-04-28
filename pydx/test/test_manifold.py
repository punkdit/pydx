
from random import random

from pydx.mjet import MJet
from pydx.manifold import RManifold
from pydx.tensor import Tensor
from pydx.field import TensorField
from pydx.test.test_field import LinearTransform, MobiusTransform


def test_rmanifold():
    n_evals = 0

    id_dd = TensorField.identity( (Tensor.dn,Tensor.dn), 2 )
    id_uu = TensorField.identity( (Tensor.up,Tensor.up), 2 )

#    T = MobiusTransform(*[((1.2+i),(3.3*i)) for i in (-2,3,1,7)])
    T = LinearTransform(*[random()-0.5 for _ in range(4)])
    g = id_dd.transform(T)
    g_uu = id_uu.transform(T)

    manifold = RManifold(g, g_uu)
#    print manifold.kretschman(0,0)[()][()]
#    print manifold.curvature(0,0)[()][()]

    assert manifold.gamma.valence == (Tensor.up,Tensor.dn,Tensor.dn)

    def test_equ(a, b, n=10):
        assert a.valence == b.valence
        xys = [(random()-0.5,random()-0.5) for _ in range(n)]
        n_evals = 0
        for x,y in xys:
            ta = a(x,y)
            tb = b(x,y)
            for idxs in ta.genidx():
                r1 = ta[idxs]
                r1 = MJet(0).promote(r1)[()]
                r2 = tb[idxs]
                r2 = MJet(0).promote(r2)[()]
                assert abs(r1-r2)<1e-10, (r1,r2)
                n_evals += 1
        return n_evals

    t_id = Tensor.identity( (Tensor.up,Tensor.dn), 2 )

    # XX This is in test_tensor
    xys = [(random()-0.5,random()-0.5) for _ in range(10)]
    for x,y in xys:
        p1 = T.partial(x,y)
        Tx,Ty=T(x,y)
        p2 = T.inverse.partial(Tx,Ty)
        p = p1.mul( p2, (1,0) )
        assert p.is_close(t_id)
        n_evals += 1

    t_id = Tensor.identity( (Tensor.up,Tensor.up), 2 )

    for x,y in xys:
        assert t_id.is_close( id_uu(x,y) )
        n_evals += 1

    g1 = g.uu.transform(T.inverse)
    n_evals += test_equ( g1, id_uu )

    n_evals += test_equ(
        TensorField.identity( (Tensor.up,Tensor.dn), 2 ),
        g.uu.mul( g.dd, (0,0) ))

    def test_sym(tfield, sym, anti=1.0, n=10):
        xys = [(random()-0.5,random()-0.5) for _ in range(n)]
        n_evals = 0
        for x,y in xys:
            tensor = tfield(x,y)
            tensor_t = tensor.transpose(*sym)
            for idxs in tensor.genidx():
                r1 = tensor[idxs][()]
                r2 = anti * tensor_t[idxs][()]
                assert abs(r1-r2)<1e-10, (r1,r2)
                n_evals += 1
        return n_evals
    # True for any manifold
    n_evals += test_sym( g.uu, (1,0) )
    n_evals += test_sym( g.dd, (1,0) )
    n_evals += test_sym( manifold.gamma, (0,2,1) )
    n_evals += test_sym( manifold.riemann, (0,1,3,2), -1.0 )
    # True in flat space
    n_evals = test_equ( manifold.riemann, TensorField.zero( manifold.riemann.valence, 2 ))

    # test for zero scalar curvature and zero kretschman
    for x, y in xys:
        assert abs(manifold.curvature(x,y)[()][()]) < 1e-10
        assert abs(manifold.kretschman(x,y)[()][()]) < 1e-10
        n_evals += 1

    print "test_rmanifold: n_evals:", n_evals

