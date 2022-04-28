#!/usr/bin/env python

#    __init__.py
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


import sys


###############################################################################

class Scalar(object):
    stack = []
#    clients = [mjet, MJet, Tensor, TensorField] # also: RManifold, Complex, ODE
    clients = [] # also: RManifold, Complex, ODE
    current = None
    def __init__( self, type, promote, zero, one, callback=None, cb_args=(), restore_cb=None ):
        self.type = type
        self.promote = promote
        self.zero = zero
        self.one = one
        self.callback = callback
        self.cb_args = cb_args 
        self.restore_cb = restore_cb

    def notify_client(self, client):
        client.scalar = self
        # deprecated:
        client.scalar_type = self.type
        client.scalar_promote = self.promote
        client.scalar_zero = self.zero
        client.scalar_one = self.one
    
    def set( cls, scalar_type, scalar_promote, scalar_zero, scalar_one, callback=None, cb_args=(), restore_cb=None ):
        cls.stack.append( cls.current )
        scalar = Scalar(scalar_type, scalar_promote, scalar_zero, scalar_one, callback, cb_args, restore_cb)
        #print "set", scalar.type
        if scalar.callback is not None:
            scalar.callback( *scalar.cb_args )
        cls.current = scalar
        for client in cls.clients:
            scalar.notify_client(client)
    set = classmethod(set)
    
    def restore(cls):
        scalar = cls.current
        if scalar.restore_cb is not None:
            scalar.restore_cb( *scalar.cb_args )
        scalar = cls.stack.pop()
        #print "restore", scalar.type, scalar.restore_cb
        if scalar.callback is not None:
            scalar.callback( *scalar.cb_args )
        for client in cls.clients:
            client.scalar_type = scalar.type
            client.scalar_promote = scalar.promote
            client.scalar_zero = scalar.zero
            client.scalar_one = scalar.one
        cls.current = scalar
    restore = classmethod(restore)

    def add_client(cls, client):
        if client not in cls.clients:
            cls.clients.append(client)
            if cls.current is not None:
                cls.current.notify_client(client)
    add_client = classmethod(add_client)

Scalar.current = Scalar(float, float, 0., 1.) # XX clients miss this... hook into clients.append

def set_scalar( *args, **kw ):
    Scalar.set(*args, **kw)
def restore_scalar():
    Scalar.restore()

def set_mpf_scalar(prec=53):
    import gmpy
    scalar_zero = gmpy.mpf(0.0)
    scalar_one = gmpy.mpf(1.0)
    scalar_type = type(scalar_zero)
    scalar_promote = gmpy.mpf
    set_scalar( scalar_type, scalar_promote, scalar_zero, scalar_one,
        gmpy.set_minprec, (prec,) )

def set_interval_scalar(prec=53):
    from mpfi import Interval, set_default_prec
    scalar_zero = Interval(0.0)
    scalar_one = Interval(1.0)
    scalar_type = type(scalar_zero)
    scalar_promote = Interval
    set_scalar( scalar_type, scalar_promote, scalar_zero, scalar_one,
        set_default_prec, (prec,) )

def set_symbolic_scalar():
    from symbolic import Float, ZeroFloat, OneFloat
    # push current scalar onto Float's "scalar" stack 
    # See eg. compile method                          
    Float.push_scalar( 
        Scalar.current.zero, 
        Scalar.current.one, 
        Scalar.current.promote, 
        Scalar.current.type, 
    )
    scalar_type = Float
    scalar_promote = Float.promote
    scalar_zero = ZeroFloat()
    scalar_one = OneFloat()
    set_scalar( scalar_type, scalar_promote, scalar_zero, scalar_one, 
        restore_cb = Float.restore_scalar )

