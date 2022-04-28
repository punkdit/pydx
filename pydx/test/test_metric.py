
from pydx.mjet import MJet
from pydx.scalar import Scalar, set_interval_scalar, restore_scalar
from pydx.scalar.mpfi import Interval
from pydx.options import options
from pydx.metric import HyperbolicMetric, SwartzchildMetric, RZCurzonMetric
from pydx.manifold import RManifold
from pydx.geodesic import Geodesic

def test_hyperbolic_geodesic():
    set_interval_scalar()

    scalar_promote = MJet.scalar_promote
    scalar_one = MJet.scalar_one
    scalrand = lambda : scalar_promote(random())

    n_evals = 0

    # Test these geodesics are circles centred
    # on a point on the x-axis.
    g = HyperbolicMetric()
    manifold = RManifold( g )
            
    x, y = scalar_promote(0.0), scalar_promote(1.0)
    dx, dy = scalar_promote(1.0), scalar_promote(0.0)
    t0 = scalar_promote(0.0)

    ode = Geodesic( manifold, t0, [x,y,dx,dy] )
    h = scalar_promote(0.01)
    for i in range(20):
        x = ode.x
        Y0 = ode.contract1( x.hull(x+h), ode.y )
        ode.istep2( h, Y0 )
        x, y, dx, dy = ode.y
        print (x**2+y**2)
        assert (x**2+y**2).contains( Interval(1.0) )
        n_evals += 1

    restore_scalar()
    print "test_hyperbolic_geodesic: n_evals:", n_evals

def test_swartzchild_geodesic():
    set_interval_scalar(options.prec)

    scalar_promote = MJet.scalar_promote
    scalar_one = MJet.scalar_one
    scalrand = lambda : scalar_promote(random())

    n_evals = 0

    M = 1.0 # mass
    g = SwartzchildMetric(M)
    manifold = RManifold( g )

    t = Interval(0.0)
    a = 1.0 # arbitrary
    dt = Interval(a)
    r = Interval(3.0)
    dr = Interval(0.0)
    theta = Interval(0.0).get_pi() / 2 # blows at sin(theta)=0
    dtheta = Interval(0.0)
    phi = Interval(0.0)
    dphi = a / ( M * Interval(27.0).sqrt() )
    x = [t, r, theta, phi, dt, dr, dtheta, dphi]
            
    t0 = scalar_promote(0.0)

    ode = Geodesic( manifold, t0, x )
#    h = scalar_promote(0.2)
    h = scalar_promote(options.step_size)
    for i in range(options.steps):
#    while phi.lower < (3.1416*2):
        x = ode.x
        Y0 = ode.contract1( x.hull(x+h), ode.y )
        ode.istepn( h, Y0, options.order )

        t, r, theta, phi, dt, dr, dtheta, dphi = ode.y
        w = sum( yi.width() for yi in ode.y )
        print t, r, phi, "width=", w
        assert r.contains( 3.0 )
        n_evals += 1

    restore_scalar()
    print "test_swartzchild_geodesic: n_evals:", n_evals

def test_rzcurzon_geodesic():
    set_interval_scalar(options.prec)

    scalar_promote = MJet.scalar_promote
    scalar_one = MJet.scalar_one
    scalrand = lambda : scalar_promote(random())

    n_evals = 0

    m = 1.0
    g = RZCurzonMetric(m)
    manifold = RManifold( g )
            
    # set from options, or set from default starting values
    r = options.get('r', -0.01)
    z = options.get('z', 0.1)
    dr = options.get('dr', -1.0)
    dz = options.get('dz', -0.2)
    t0 = options.get('t0', 0.0)
    h = options.get('step_size', 0.0001)
    r,z,dr,dz,t0,h = [scalar_promote(_) for _ in r,z,dr,dz,t0,h ]
    # XX should we just keep all of those as params in options ? XX

    # build the ode object:
    ode = Geodesic( manifold, t0, [r,z,dr,dz], paramname='t', varnames='r z dr dz'.split() )

    file = options._file
    if options.filename is not None and options.resume:
        # replay file
        info = Info(begin=False)
        try:
            while 1:
                ode.loadstate(file)
        except:
            pass
        print "done loadstate"
        file.close()
        file = open(options.filename, 'a')
    else:
        ode.dumpstate(file)

    try:
        while ode.i<options.steps or options.steps==0:
            x = ode.x
            Y0 = ode.contract1( x.hull(x+h), ode.y ) #, max_iter=5 )
            # ode.istep2( h, Y0 )
            ode.istepn( h, Y0, options.order )
            ode.dumpstate(file)
            width = sum( yi.width() for yi in ode.y )
            abs_size = sum( abs(yi).upper for yi in ode.y )
            abs_size = max( abs_size, 1e-12 )
            relwidth = width / abs_size
            print >>file, 'width =', repr(width)
            print "step", ode.i
            if relwidth > 1e5:
                print >>file, 'relwidth=%r'%relwidth
                break
    except DisorderException, e:
        print >>file, 'e="%s"' % repr(e)
    print "Finished"

    restore_scalar()

#
#
# XX another test to write: that as we increase the order,   XX
# XX solutions overlap each other                            XX
#

