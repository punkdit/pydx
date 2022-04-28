#!/usr/bin/env python

#    optimize.py
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
import os
from os.path import join as pjoin
from time import time

import gmpy
from gmpy import mpf
mpf_type = type(mpf(1.0))

import mpfi
from mpfi import Interval as _Interval
from mpfi import get_default_prec, set_default_prec
from mpfi import promote as _promote
from mpfi import from_string as interval_from_string

import ctypes

from pyjit import llvm
from pyjit import jit

filename = pjoin(os.getcwd(),mpfi.__file__) 
jit.load_lib( filename )

print filename
#open(filename)
#mpfi_lib = ctypes.cdll.LibraryLoader( filename )
#mpfi_lib = ctypes.cdll.LibraryLoader( mpfi.__file__ )

import dl
mpfi_lib = dl.open( filename )


def lldecl():
    decls = []
    f = open("_mpfi.h")
    for line in f.readlines():
        line = line.strip()
        assert line.startswith('int ')
        assert line.endswith(');')
        line = line[ len('int') : -len(';') ]
        line = line.strip()
        func, rest = line.split('(')
        func = func.strip()
        n_args = rest.count(',')+1
        args = [ llvm.ubyte.pointer ] * n_args
        decl = llvm.Declare( func, llvm.int_, args )
        decls.append( decl )
#    # mpfi_pow
#    bb = llvm.BasicBlock()
#    a_v=llvm.Var()
#    b_v=llvm.Var()
#    res_v=llvm.Var()
#    tmp=llvm.Var()
#    bb.append( llvm.Call( llvm.int_, "mpfi_log", [(llvm.ubyte.pointer,tmp_v),(llvm.ubyte.pointer,a_v)]))
#    bb.append( llvm.Call( llvm.int_, "mpfi_mul",
#        [(llvm.ubyte.pointer,res_v),(llvm.ubyte.pointer,arg0_v),(llvm.ubyte.pointer,arg1_v)]))
#    f = llvm.Function( "mpfi_pow", llvm.void,
#        [ (llvm.ubyte.pointer,res_v), (llvm.ubyte.pointer,a_v), (llvm.ubyte.pointer,b_v) ], [bb] )
#    print f
#    blow
#    decls.append( f )
    return decls

def init_lldecl():
    global __mc_decl
    decls = lldecl()
    m = llvm.Module( decls )
#    print m.llstr()
    __mc_decl = jit.consume( m.llstr(), verify=True, dump=True )

init_lldecl()

def split( x, n=2 ):
    if len(x)==0:
        return [x]
    if len(x)==1:
        return [ [_x] for _x in x[0].split(n) ]
    x2s = split( x[1:], n )
    return [ [x1]+x2 for x2 in x2s for x1 in x[0].split(n) ]

#print split( (_Interval(0,1),) )
#print split( (_Interval(0,1),_Interval(0,1)) )
#sys.exit()

#def minimize( func, freevars, EPSILON=1e-2, n=10 ):
#    """ 
#        func: python function object, vs: list of Interval objects
#        minimize func on domain given by vs
#        return x, y, x:list of no more than about n intervals, y: interval ( minimum )
#    """
#
##    try:
##        func = eval(code)
##    except MemoryError:
##        print code
##        raise
#
#    if len(freevars)==0:
#        # no free vars
#        result = _Interval(func({}))
#        return [], result
#
#    vs = freevars.values()
#    print len(vs)
#    labels = [v.label for v in vs]
#    xss = [ [v.x for v in vs] ] # list of args to func
#    while 1:
##        for x in xs:
##            print "minimize  :", ','.join([str(xc) for xc in x]),'=', func(*x)
#
#        ps = [ (xs,func(xs)) for xs in xss ]
#        ps.sort( lambda p, q: cmp(p[1].lower,q[1].lower) )
#        idx = 1
#        while idx<len(ps):
#            if ps[0][1].upper < ps[idx][1].lower:
#                ps.pop(idx)
#            else:
#                idx += 1
#        w = min( [ p[1].width() for p in ps ] )
#        if w<EPSILON or len(ps)>n:
#            break
#        print "minimize",len(xss), w
#
#        _ps = ps[:]
#        _ps.sort( lambda a,b : cmp( b[1].width(), a[1].width() ) )
#        p = _ps.pop(0)
#
#        xss = []
#        for xs,_ in ps:
#            _xss = split(xs,2)
#            xss += _xss
#        print "minimize",len(xss), w
#        print "."
#    xss = [p[0] for p in ps]
#    ys = [p[1] for p in ps]
#    lower = ys[0].lower
#    upper = min( [ y.upper for y in ys ] )
##    assert type(lower)==mpf_type
##    assert type(upper)==mpf_type
#    result = _Interval(lower,upper)
#    return xss, result

def split_greatest_var( func, xs ):
    x0s = [ _Interval(x.centre()) for x in xs ]
    ys = []
    for i in range(len(xs)):
        _x = x0s[i]
        x0s[i] = xs[i]
        ys.append( func( x0s ) )
    idx = ys.index( max(ys) )
    left,right = xs[idx].bisect()
    x1s = xs[:]
    x2s = xs[:]
    x1s[idx] = left
    x2s[idx] = right
    return x1s, x2s

def minimize( func, freevars, EPSILON=1e-2, n=500 ):
    """ 
        func: python function object, vs: list of Interval objects
        minimize func on domain given by vs
        return x, y, x:list of no more than about n intervals, y: interval ( minimum )
    """

    if len(freevars)==0:
        # no free vars
        result = _Interval(func([]))
        return [], result

    vs = freevars.values()
    vs.sort( lambda a,b: cmp(a.label,b.label) )
#    print len(vs)
    labels = [v.label for v in vs]
    xss = [ [v.x for v in vs] ] # list of args to func, a list of "boxes"
    while 1:
        ps = [ (xs,func(xs)) for xs in xss ]
        ps.sort( lambda p, q: cmp(p[1].lower,q[1].lower) )
        idx = 1
        while idx<len(ps):
            if ps[0][1].upper < ps[idx][1].lower:
                ps.pop(idx)
            else:
                idx += 1
        w = min( 
            [ min( p[1].width()/max(1e-6,abs(p[1].centre())), p[1].width() )
                for p in ps ] )
        if w<EPSILON or len(ps)>n:
            break
#        print "minimize",len(ps), w

        xss = []
        for xs,_ in ps:
            x1s, x2s = split_greatest_var(func, xs)
            xss.append( x1s )
            xss.append( x2s )
#        print "minimize",len(xss), w
#    print "------->",len(ps), w

    xss = [p[0] for p in ps]
    ys = [p[1] for p in ps]
    lower = ys[0].lower
    upper = min( [ y.upper for y in ys ] )
#    assert type(lower)==mpf_type
#    assert type(upper)==mpf_type
    result = _Interval(lower,upper)
    return xss, result

#
###############################################################################
#

class Interval(object):
    uid = 0
    file = open('dump.ll','w')
    def __init__( self, lower=None, upper=None ):
        self.freevars = {} # map label (str) to interval
        # cache lower and upper calculations:
        self._lower = None
        self._upper = None
        self.funcname = None
        if lower is not None:
            if isinstance(lower,Interval):
                assert upper is None
                lower = lower.x
#            if type(lower) in (int,long,float) and upper is None:
#                self.label = str(lower)
#                self.x = lower
#                return
            if hasattr(lower,"__float__") or type(lower)==mpf_type:
                lower = _Interval(lower,upper)
                upper = None
            if isinstance(lower,_Interval):
                assert upper is None
                if lower.width()==0.0:
#                    self.label = "interval_from_string('%s')"%lower.to_string()
#                    self.label = str(lower.centre())
                    self.label = lower.centre()
                    self.x = lower # we still set this
                else:
                    self.x = lower # domain of this variable
                    self.label = "v%s"%Interval.uid
                    Interval.uid += 1
                    self.freevars[self.label] = self
            else:
                assert not isinstance(lower,Interval), repr(lower)
                assert 0, repr(lower)
            
    def __str__( self ):
        return "[[%.7f,%.7f]]" % (self.lower, self.upper)
    def __repr__( self ):
#        return "%s(%s)"%( self.expr(), [(key,v.lower,v.upper) for key,v in self.freevars.items()] )
#        return "%s(%.16f,%.16f)" % (self.__class__.__name__,self.lower, self.upper)
        return "%s(%.16f,%.16f)" % (self.__class__,self.lower, self.upper)

    def freeze( self ):
        return Interval( self.lower, self.upper )

    def expr( self ):
        return repr(self.label)

    def promote(cls,x):
        if isinstance(x,cls):
            return x
        x = _promote(x)
        if x is None:
            return None
        return Interval(x)
    promote = classmethod(promote)

    def deeplen( self ):
        return 1

    def genll( self ):
        # this uses a 32 byte alloca per node, so we could blow the stack...
        const_cache = {} # map addr to Var's
        var_cache = [] # list of alloca'ed Var's no longer needed
        bb = llvm.BasicBlock()
        arg_v = self._genll( bb, const_cache, var_cache )
        res_v = llvm.Var()
        bb.append(
            llvm.Call( llvm.int_, "mpfi_set", 
                [(llvm.ubyte.pointer,res_v),(llvm.ubyte.pointer,arg_v)]))
        self.free_var( var_cache, arg_v )
        for var in var_cache:
            bb.append( llvm.Call( llvm.int_, "mpfi_clear", [(llvm.ubyte.pointer,var)]))
        bb.append( llvm.Ret() )
        args = [ (llvm.ubyte.pointer,res_v) ]
        labels = self.freevars.keys()
        labels.sort()
        args += [ (llvm.ubyte.pointer,label) for label in labels ]
        f = llvm.Function(None, llvm.void, args, [bb])
        return f

    def alloc_var( self, bb, var_cache ):
        if var_cache:
            return var_cache.pop(0)
        res_v = llvm.Var()
        bb.append( llvm.Assign( res_v, llvm.Alloca( llvm.ubyte, mpfi.sizeof_mpfi ) ))
        bb.append( llvm.Call( llvm.int_, "mpfi_init", [(llvm.ubyte.pointer,res_v)]))
        return res_v

    def free_var( self, var_cache, var ):
        assert isinstance(var,llvm.Var)
        var_cache.append( var )

    def _genll( self, bb, const_cache, var_cache, *args ):
        " this is overloaded for each node type. return var that we assigned to. "
        if type(self.label)==mpf_type:
            addr = self.x.get_addr()
#            if addr in const_cache:
#                return const_cache[addr]
            addr_v = llvm.Var()
            tmp_v = llvm.Var()
            bb.append( llvm.Assign( tmp_v, llvm.Cast( llvm.long_, addr, llvm.ubyte.pointer ) ))
        else:
            tmp_v = llvm.Var.promote( self.label )
        res_v = self.alloc_var( bb, var_cache )
        bb.append( llvm.Call( llvm.int_, "mpfi_set", [(llvm.ubyte.pointer,res_v),(llvm.ubyte.pointer,tmp_v)]))
#        if type(self.label)==mpf_type:
#            const_cache[ addr ] = res_v
        return res_v

    def get_iop( self, func_args ):
#        print "argnum = ", func_args.index(self), self, func_args
        for argnum,arg in enumerate(func_args):
            if arg is self:
                iop = mpfi.iop_new( argnum = argnum )
                return iop
        assert 0

    def __call__( self, xs ):
        if self.funcname is None:
            f = self.genll()
#            print self.freevars.keys()
#            print f
#            print "building..."
#            sys.stdout.flush()
            s = f.llstr()
            if self.file:
                print >>self.file, s
                self.file.flush()
#            t0 = time()
#            print f.name, len(f.bbs[0])
#            print f
            jit.consume( s, verify=True, dump=True ) #, f.name )
            self.funcname = f.name
#            print 'time\t\t', time()-t0
#            print "done"
            
        addr = jit.get_func(self.funcname)
    
        labels = self.freevars.keys()
        labels.sort()
        args = ( ctypes.c_void_p, ) * ( len( labels ) + 1 )
        func_ptr_type = ctypes.CFUNCTYPE(None, *args)
        
        func_ptr = func_ptr_type(addr)
        
        res = _Interval()
        args = [ ctypes.c_void_p( res.get_addr() ) ]
        args += [ ctypes.c_void_p( xs[i].get_addr() ) for i in range(len(xs)) ]
        func_ptr( *args )
        return res

    def __neg__( x ):
        return Neg(x)
    def __pos__( x ):
        return x
#        return Pos(x)
    def __abs__( x ):
        return Abs(x)
    def __add__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Add(x,x.cast(y))
    def __sub__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Sub(x,x.cast(y))
    def __mul__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Mul(x,x.cast(y))
    def __div__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Div(x,x.cast(y))
    def __radd__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Add(x.cast(y),x)
    def __rsub__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Sub(x.cast(y),x)
    def __rmul__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Mul(x.cast(y),x)
    def __rdiv__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Div(x.cast(y),x)
    def __pow__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return NotImplemented
        return Exp(Mul(x.cast(y),Log(x)))
        return Pow(x,x.cast(y))

    def cast( self, x ):
        if isinstance(x, Interval):
            return x
        if not isinstance(x, _Interval):
            x = _Interval(x)
        return Interval(x)

#    def __coerce__( self, other ):
#        print "Interval.__coerce__", self, other
#        if hasattr(other,"__float__"):
#            return (self, Const(other))
#        other = promote(other)
#        if other is not None:
#            return (self,Interval(other))
#        return None

#    def get_vars( self ):
#        vs = self.freevars.values()
#        # this is the arg order that we use in the call to our func:
#        vs.sort( lambda a,b: cmp(a.label,b.label) )
#        return vs

    def get_lower( self ):
#        print "\nminimize( %s )" % (self.code())
        if self._lower is not None:
            lower = self._lower
        else:
            xs, y = minimize( self, self.freevars )
            lower = y.lower
            self._lower = lower
#        print "minimize lower =", lower
        return lower
    lower = property(get_lower)
    def get_upper( self ):
#        print "\nminimize( %s )" % ((-self).code())
        if self._upper is not None:
            upper = self._upper
        else:
            xs, y = minimize( -self, self.freevars )
            upper = -y.lower
            self._upper = upper
#        print "minimize upper =", upper
        assert upper >= self.lower, (upper, self.lower)
        return upper
    upper = property(get_upper)


    def width(x):
        return x.upper - x.lower
    def overlapping(x, y):
        return x.upper >= y.lower

    def intersect(x, y):
        return Interval(max(x.lower, y.lower), min(x.upper, y.upper))
    def hull(x, y):
        return Interval(min(x.lower, y.lower), max(x.upper, y.upper))
    def contains(self, other):
        other = self.cast(other)
        return self.lower <= other.lower and self.upper >= other.upper
    def inside(self, other):
        other = self.cast(other)
        return other.lower <= self.lower and other.upper >= self.upper
    def centre(self):
        return (self.lower+self.upper)/2.0

    def __gt__( x, y ):
        y = Interval.promote(y)
        if x.lower > y.upper:
            return True
        return False
    def __ge__( x, y ):
        y = Interval.promote(y)
        if x.lower > y.lower:
            return True
        if x == y:
            return True
        return False
    def __lt__( x, y ):
        y = Interval.promote(y)
        if x.upper < y.lower:
            return True
        return False
    def __le__( x, y ):
        y = Interval.promote(y)
        if x.upper < y.lower:
            return True
        if x == y:
            return True
        return False
    def __ne__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return True
        if x.lower != y.lower or x.upper != y.upper:
            return True
        return False
    def __eq__( x, y ):
        y = Interval.promote(y)
        if y is None:
            return False
        if x.lower == y.lower and x.upper == y.upper:
            return True
        return False


#    def is_close(x, y):
#    def bisect(x):
    
#    def ceil(x):
#    def floor(x):
    
#    def max(x, y):
#    def min(x, y):
    
#    def sqr(x):
    
    #// transcendental increasing monotonic
    
#    def acosh(x):
#    def asin(x):
#    def asinh(x):
#    def atan(x):
#    def atanh(x):
    def exp(x):
        return Exp(x)
#    def pow(x, y):
#    def log(x):
#    def sinh(x):
#    def sqrt(x):
#    def tanh(x):
    
    #// transcendental decreasing monotonic
    
#    def acos(x):
    
    #// transcendental general
    
#    def cosh(x):
#    
#    def cos(x):
#    
#    def sin(x):
#
#    def tan(x):

#class Const(Interval):
#    def __init__( self, x ):
#        assert isinstance(x,_Interval)
#        Interval.__init__(self)
#        self.x = x
#    def expr( self ):
##        print repr(self.x)
#        return "interval_from_string('%s')"%self.x.to_string() # _Interval(...)
#    def str( self ):
#        return str(self.x)

#def promote(value):
#    if hasattr(value,"__float__"):
#        return Const(_Interval(value))
#    if isinstance(value,_Interval):
#        if value.width()==0.0:
#            return Const(value)
#        return Var(value)
#    if isinstance(value,Interval):
#        return value
#    raise TypeError, str(value)
#    # and what about Const ??

class Multi(Interval):
    def __init__( self, args ):
        Interval.__init__(self)
        for arg in args:
            assert isinstance(arg,Interval)
        self.args = args
        for arg in args:
            self.freevars.update(arg.freevars)
    def deeplen( self ):
        return 1 + sum([ arg.deeplen() for arg in self.args ])
    def get_iop( self, func_args ):
        iops = [ arg.get_iop(func_args) for arg in self.args ]
#        op_addr = ctypes.c_long()
#        op_addr.value = self.mpfi_func 
#        op_addr = eval(str(self.mpfi_func)[-11:-1]+'L') # ctypes hack to cast pointer to long :P
        op_addr = mpfi_lib.sym(self.mpfi_func)
        iop = mpfi.iop_new( op_addr, iops )
        return iop

class Unary(Multi):
    def __init__( self, x ):
        Multi.__init__(self,[x])
    def _genll( self, bb, funcname, const_cache, var_cache, *args ):
        res_v = self.alloc_var( bb, var_cache )
        if self.funcname:
#            print "got one!", id(self)
            args = [(llvm.ubyte.pointer,res_v)]
            labels = self.freevars.keys()
            labels.sort()
            args += [(llvm.ubyte.pointer,label) for label in labels]
            bb.append( llvm.Call( llvm.void, self.funcname, args ) )
        else:
            arg_v = self.args[0]._genll(bb,const_cache, var_cache,*args)
            bb.append( llvm.Call( llvm.int_, funcname, [(llvm.ubyte.pointer,res_v),(llvm.ubyte.pointer,arg_v)]))
            self.free_var( var_cache, arg_v )
        return res_v

class Binary(Multi):
    def __init__(self,x,y):
        Multi.__init__(self,[x,y])
    def _genll( self, bb, funcname, const_cache, var_cache, *args ):
#        print repr(self.args)
        res_v = self.alloc_var( bb, var_cache )
        if self.funcname:
#            print "got one!", id(self)
            args = [(llvm.ubyte.pointer,res_v)]
            labels = self.freevars.keys()
            labels.sort()
            args += [(llvm.ubyte.pointer,label) for label in labels]
            bb.append( llvm.Call( llvm.void, self.funcname, args ) )
        else:
            arg0_v = self.args[0]._genll(bb,const_cache, var_cache,*args)
            arg1_v = self.args[1]._genll(bb,const_cache, var_cache,*args)
    #    res_v = llvm.Var()
    #    bb.append( llvm.Assign( res_v, llvm.Alloca( llvm.ubyte, mpfi.sizeof_mpfi ) ))
    #    bb.append( llvm.Call( llvm.int_, "mpfi_init", [(llvm.ubyte.pointer,res_v)]))
            bb.append( llvm.Call( llvm.int_, funcname,
                [(llvm.ubyte.pointer,res_v),(llvm.ubyte.pointer,arg0_v),(llvm.ubyte.pointer,arg1_v)]))
            self.free_var( var_cache, arg0_v )
            self.free_var( var_cache, arg1_v )
        return res_v

class Neg(Unary):
    mpfi_func = "mpfi_neg"
    def expr( self ):
        return "(-%s)"%self.args[0].expr()
    def str( self ):
        return "(-%s)"%self.args[0].str()
    def _genll( self, bb, *args ):
        return Unary._genll( self, bb, "mpfi_neg", *args )
#    def __call__( self, kw={} ):
#        return -self.args[0](kw)
#class Pos(Unary):
#    def expr( self ):
#        return "(+%s)"%self.args[0].expr()
#    def str( self ):
#        return "(+%s)"%self.args[0].str()
#    def _genll( self, bb, *args ):
#        arg_v = self.args[0]._genll(bb, *args)
#        return arg_v
##    def __call__( self, kw={} ):
##        return self.args[0](kw)
class Abs(Unary):
    mpfi_func = "mpfi_abs"
    def expr( self ):
#        return "(%s.abs())"%self.args[0].expr()
        return "abs(%s)"%self.args[0].expr()
    def str( self ):
        return "abs(%s)"%self.args[0].str()
    def _genll( self, bb, *args ):
        return Unary._genll( self, bb, "mpfi_abs", *args )
#    def __call__( self, kw={} ):
#        return abs(self.args[0](kw))
class Exp(Unary):
    mpfi_func = "mpfi_exp"
    def expr( self ):
        return "(%s.exp())"%self.args[0].expr()
    def str( self ):
        return "exp(%s)"%self.args[0].str()
    def _genll( self, bb, *args ):
        return Unary._genll( self, bb, "mpfi_exp", *args )
#    def __call__( self, kw={} ):
#        return (self.args[0](kw)).exp()
class Log(Unary):
    mpfi_func = "mpfi_log"
    def expr( self ):
        return "(%s.log())"%self.args[0].expr()
    def str( self ):
        return "log(%s)"%self.args[0].str()
    def _genll( self, bb, *args ):
        return Unary._genll( self, bb, "mpfi_log", *args )


class Add(Binary):
    mpfi_func = "mpfi_add"
    def expr( self ):
        return "(%s+%s)"%tuple([arg.expr() for arg in self.args])
    def str( self ):
        return "(%s+%s)"%tuple([arg.str() for arg in self.args])
    def _genll( self, bb, *args ):
        return Binary._genll( self, bb, "mpfi_add", *args )
#    def __call__( self, kw={} ):
#        return self.args[0](kw) + self.args[1](kw)
class Sub(Binary):
    mpfi_func = "mpfi_sub"
    def expr( self ):
        return "(%s-%s)"%tuple([arg.expr() for arg in self.args])
    def str( self ):
        return "(%s-%s)"%tuple([arg.str() for arg in self.args])
    def _genll( self, bb, *args ):
        return Binary._genll( self, bb, "mpfi_sub", *args )
#    def __call__( self, kw={} ):
#        return self.args[0](kw) - self.args[1](kw)
class Mul(Binary):
    mpfi_func = "mpfi_mul"
    def expr( self ):
        return "(%s*%s)"%tuple([arg.expr() for arg in self.args])
    def str( self ):
        return "(%s*%s)"%tuple([arg.str() for arg in self.args])
    def _genll( self, bb, *args ):
        return Binary._genll( self, bb, "mpfi_mul", *args )
#    def __call__( self, kw={} ):
#        return self.args[0](kw) * self.args[1](kw)
class Div(Binary):
    mpfi_func = "mpfi_div"
    def expr( self ):
        return "(%s/%s)"%tuple([arg.expr() for arg in self.args])
    def str( self ):
        return "(%s/%s)"%tuple([arg.str() for arg in self.args])
    def _genll( self, bb, *args ):
        return Binary._genll( self, bb, "mpfi_div", *args )
#    def __call__( self, kw={} ):
#        return self.args[0](kw) / self.args[1](kw)
class Pow(Binary):
    def expr( self ):
        return "(%s**%s)"%tuple([arg.expr() for arg in self.args])
    def str( self ):
        return "(%s**%s)"%tuple([arg.str() for arg in self.args])
    def _genll( self, bb, *args ):
        assert 0, (self.args)
        # use mpfi_exp ?
#    def __call__( self, kw={} ):
#        return self.args[0](kw) ** self.args[1](kw)

#
###############################################################################
#

def test_00():
    x = _Interval(0.0,1.0)
    _y = 3.0 * x*(1-x) + 4.0 + _Interval(0.0,0.1)
    x = Interval(x)
    y = 3.0 * x*(1-x) + 4.0 + _Interval(0.0,0.1)
    print _y, ">>", y
    assert y.inside(_y)

    x1 = _Interval(0,1)
    x2 = _Interval(0,1)
    _y = 3.0 * (x1-x2) * x1 + (x2-x1)*x2 + x2*_Interval(0,1)
    print 
    _x1 = Interval(x1)
    _x2 = Interval(x2)
    y = 3.0 * (_x1-_x2) * _x1 + (_x2-_x1)*_x2 + x2*_Interval(0,1)
    print _y, ">>", y
    assert y.inside(_y)
    print
    return y

def test_01():
    x1 = _Interval(0,1)
    x2 = _Interval(0,1)
    x3 = _Interval()
    print x1+x2

    bb = llvm.BasicBlock()
    bb.append(
        llvm.Assign('x', 
            llvm.Call( llvm.int_, "mpfi_add", 
                [(llvm.ubyte.pointer,'x3'), (llvm.ubyte.pointer,'x1'), (llvm.ubyte.pointer,'x2'), ] ) ))
    call = llvm.Call( llvm.int_, "mpfi_add", [(llvm.ubyte.pointer,'x3'), (llvm.ubyte.pointer,'x1'), (llvm.ubyte.pointer,'x2'), ] )
    print call.ty
    print "%s* "%call.ty and call.ty or ""

    bb.append( llvm.Ret(llvm.int_,'x') )
    f = llvm.Function("add", llvm.int_, [(llvm.ubyte.pointer,'x3'), (llvm.ubyte.pointer,'x1'), (llvm.ubyte.pointer,'x2'), ], [bb])
    f.to_ssa()
    decl = llvm.Declare( "mpfi_add", llvm.int_, (llvm.ubyte.pointer,llvm.ubyte.pointer,llvm.ubyte.pointer) )
    m = llvm.Module( (decl,f) )
    print m

    from pyjit import jit
    jit.init_jitcore()
    jit.load_lib( mpfi.__file__ )

    mc = jit.consume( m.llstr() )
    addr = mc.get_func(f.name)
    print "addr", addr

    import ctypes

    # typedef int (*func_ptr_type)(char *s);
    func_ptr_type = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p )
    
    # func_ptr = (func_ptr_type)123;
    func_ptr = func_ptr_type(addr)
    
    print "x3=", x3
    # x = func_ptr('arg')
    x = func_ptr( ctypes.c_void_p(x3.get_addr()), ctypes.c_void_p(x1.get_addr()), ctypes.c_void_p(x2.get_addr()), )

    print "x3=", x3

def test_02():
    # test for mem-leak
    i = 0
    while 1:
        print
        print i
        _x1 = Interval(0,1)
        _x2 = Interval(0,1)
        y = 3.0 * (_x1-_x2) * _x1 + (_x2-_x1)*_x2 + _x2*_x2

        for i in range(2):
            y = 3.0 * (y-_x2) * _x1 + (_x2-_x1)*y + _x2*_x2
        print y
        i += 1

def test_03():
    r = Interval(3.58)
    x = Interval(0.5)
    unit = Interval(1.0)
    y = x
    for i in range(10):
        y = r*y*(unit-y)
        if i%4==3:
            z=y( [r.x,x.x] )
    # r*x*(1-x), r=3.58
    print y.deeplen()
    mpfi.iop_arena_init( y.deeplen() )
    iop = y.get_iop( [r,x,unit] )
#    print iop
    sz = y.deeplen()
    TRIALS=1000
    t0 = time()
    for i in range(TRIALS):
        z,count=mpfi.iop_eval( sz, iop, [r.x,x.x,unit.x] )
    t1 = time()-t0
    print "count", count
    print "iop time:", t1/TRIALS
    print "iop_eval", z
    mpfi.iop_deepfree(iop)
    mpfi.iop_arena_free()
    print "ll:"
    t0 = time()
    for i in range(TRIALS):
        z=y( [r.x,x.x] )
    t2 = time()-t0
    print "jit time", t2/TRIALS
    print "speedup", t1/t2
    print "ll eval", z

if __name__ == "__main__":
    test_03()



