
from pydx.ode import ODE
from pydx.scalar import set_interval_scalar, restore_scalar, Scalar
from pydx.scalar.mpfi import Interval
from pydx.tensor import Tensor

def test_ode():
    set_interval_scalar()

    scalar_promote = ODE.scalar_promote
    scalar_one = ODE.scalar_one
    scalrand = lambda : ODE.scalar_promote(random())

    steps = 1000
    h = ODE.scalar_promote(0.001)

    n_evals = 0

    class Pendulum(ODE):
        def __init__( self, x0, y0 ):
            """
            """
            ODE.__init__( self, x0, y0 )
            assert len(y0)==2
        def __call__( self, x, y ):
            y, dy = y
            result = Tensor( (Tensor.up,), self.dim )
            result[0] = -dy
            result[1] = y
            return result


    # test first order interval method
    # y''(x) = -y(x), y(0) = 1, y'(0) = 0
    # solution: y(x) = cos(x)
    x = ODE.scalar_zero
    y, dy = ODE.scalar_one, ODE.scalar_zero
    ode = Pendulum( x, [y,dy] )
    for i in range(steps):
        x = Interval(ode.x).lower
        Y0 = ode.contract1( Interval(x,(x+h).upper), ode.y )
        ode.istep1( h, Y0 )
        y, dy = ode.y
#        if i%100==0:
#            print "width:", y.width(), "centre:", y.centre()
        assert ode.y[0].overlapping( ode.x.cos() )
        n_evals += 1
    print "istep1 x",ode.x," width:", ode.y[0].width(), "centre:", ode.y[0].centre()

    # test second order interval method
    # y''(x) = -y(x), y(0) = 1, y'(0) = 0
    # solution: y(x) = cos(x)
    x = ODE.scalar_zero
    y, dy = ODE.scalar_one, ODE.scalar_zero
    ode = Pendulum( x, [y,dy] )
    for i in range(steps):
        x = Interval(ode.x).lower
        Y0 = ode.contract1( Interval(x,(x+h).upper), ode.y )
        ode.istep2( h, Y0 )
        y, dy = ode.y
#        if i%100==0:
#            print "width:", y.width(), "centre:", y.centre()
        assert ode.y[0].overlapping( ode.x.cos() )
        n_evals += 1
    print "istep1 x",ode.x," width:", ode.y[0].width(), "centre:", ode.y[0].centre()

    restore_scalar()
    #assert len(Scalar.stack)==0
    print "test_ode: n_evals:", n_evals

