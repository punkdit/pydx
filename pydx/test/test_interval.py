#!/usr/bin/env python

import sys
from random import *

#from interval import Interval
#from xinterval import Interval
from pydx.scalar.mpfi import Interval
#from pydx.scalar.optimize import Interval

def test_nocrash():
    x = Interval(0,1)
    print x
    y = Interval(2,3)
    print y
    print x+y
    print x*y
    print x/y

def test_interval():
    " check that interval ops encompass the equivalent float ops "
    def f(x): return x**3 + 3*x + x*x + 99

    x = 1.0
    ix = Interval(x)
    assert f(ix).contains( f(x) ) 

    ops = [
        lambda a,b:a+b,
        lambda a,b:a-b,
        lambda a,b:a*b,
    ]

    from gmpy import mpf
    seed(1)
    for _ in range(100):
        i = 0
        x = [12345L,654321L]
        r = 1e-10
#        x = [mpf( x[0], 1024 ), mpf( x[1], 1024 )]
        y = [Interval(float(x[0])-r,float(x[0])+r),Interval(float(x[1])-r,float(x[1])+r)]
#        x = [99*(y[0]/99),99*(y[1]/99)]
#        y = [y[0].hull(x[0]),y[1].hull(x[1])]
        while 1:
            assert y[0].contains(x[0]), "\n  i=%d, \n  y=%s, \n  x=%s" % (i, repr(y), repr(x) )
            assert y[1].contains(x[1]), "\n  i=%d, \n  y=%s, \n  x=%s" % (i, repr(y), repr(x) )
            idx = choice([0,1])
            if abs(x[0])+abs(x[1])>1e10:
                if abs(x[0])>abs(x[1]):
                    op = lambda x,y:abs(x)-abs(y)
                else:
                    op = lambda x,y:abs(y)-abs(x)
#            elif abs(x[0])<1e-5:
#                op = lambda x,y:(x+y+1.0)
            else:
                op = choice(ops)
            x[idx] = op(*x)
            y[idx] = op(*y)
            w = y[0].width()+y[1].width()
#            print ops.index(op), w
            if w>1.0:
                print "break"
                print "  w =", w
                print "  i =", i
                break
            i += 1






