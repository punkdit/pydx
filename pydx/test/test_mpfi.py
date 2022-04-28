#!/usr/bin/env python


from pydx.scalar.mpfi import Interval


def test_simple():
    x = Interval(0,1)
    x = Interval(0.0,1.0)
    
    print x
    
    x = Interval(0.1,0.1000000)
    x = x.exp()
    assert x.width()<1e-14, x.width()
    
    print x



