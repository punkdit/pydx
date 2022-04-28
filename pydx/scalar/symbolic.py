#!/usr/bin/env python

#    symbolic.py
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

from random import random, seed
from time import time

import fmath

try:
    import gmpy
    from gmpy import mpf
except ImportError:
    print " ** no gmpy found ** "

try:
    from xinterval import Interval, get_default_prec, set_default_prec, DisorderException
except ImportError, e:
    print " ** could not import xinterval ** "
    print e

__module__ = sys.modules[__name__]

class ComputationContext(object):
    def __init__( self, name, args, dumpfile=None ):
        """
            @param name: function name
            @param args: function arg names (can build lvalue's out of these)
            @param dumpfile: dump the function to this file (optional)
        """
        self.lines = []
        self.append( 'def %s(%s):' % (name, ', '.join( [str(v) for v in args] )))
        self.float_cache = {}
        self.var_cache = {}
        self.var_floats = []
        self.name = name
        self.dumpfile = dumpfile
    def append( self, line ):
        self.lines.append(line)
    def assign( self, lvalue, a_float ):
        """
            @param lvalue: a string repr of the left hand side of the assign
            @param a_float: the right hand side (abstract float) of the assign
        """
        #Generate a sequence of unique "assign" statements:
        #  float_cache: map Float -> VarFloat
        #  var_cache: map VarFloat -> Float  (an "assignment")
        #  var_floats: what order to assign the vars in
        float_cache, var_cache, var_floats\
            = a_float.uniq_form(self.float_cache, self.var_cache, self.var_floats)
        assert float_cache is self.float_cache
        assert var_cache is self.var_cache
        assert var_floats is self.var_floats
        self.var_cache[lvalue] = self.float_cache[a_float] # ASSIGNMENT
        self.var_floats.append( lvalue )
    def finalize(self):
        """
            Generate the function. (exec the source.)
        """
        for var in self.var_floats:
            self.append( "  %s = %s" % (var, self.var_cache[var]) )
        self.append( "  return" )
        src = '\n'.join(self.lines)
        if self.dumpfile is not None:
            print >>self.dumpfile, src
        exec src
        return locals()[self.name]


class Float(object):
    scalar_zero = 0.0
    scalar_one = 1.0
    scalar_promote = float
    scalar_type = float
    old_scalar = None
    def push_scalar(cls, zero, one, promote, type): # XX wrong order
        assert cls.old_scalar is None, "stack depth of 1!"
        cls.old_scalar = cls.scalar_zero, cls.scalar_one, cls.scalar_promote, cls.scalar_type # XX wrong order
        cls.scalar_zero = zero
        cls.scalar_one = one
        cls.scalar_promote = promote
        cls.scalar_type = type
        #print "Float.push_scalar", cls.scalar_type
    push_scalar = classmethod(push_scalar)
    def restore_scalar(cls):
        cls.scalar_zero, cls.scalar_one, cls.scalar_promote, cls.scalar_type = cls.old_scalar # XX wrong order
        #print "Float.restore_scalar", cls.scalar_type
        cls.old_scalar = None
    restore_scalar = classmethod(restore_scalar)

    cache = {} # Fly-weight design pattern
    unop_cache = {} # map name to class, eg. SinFloat

    def __new__(cls, val=None):
        ob = object.__new__(cls)
        ob._hash = None
        assert isinstance(val,cls.scalar_type), val
        ob.val = val
        if ob.val == cls.scalar_zero:
            assert ZeroFloat._zero_float is None # singleton
        if ob in cls.cache:
            return cls.cache[ob]
        cls.cache[ob] = ob
        return ob
#    def __init__(self, val=None):
#        pass
    def __init__(self, *args):
        pass
    def __hash__(self):
        if self._hash is None:
            self._hash = self.hash()
        #assert type(self._hash)==int, self._hash
        return self._hash
    def hash(self):
        return hash(self.val)
    def __eq__(self, other):
        return type(self)==type(other) and self.val==other.val
    def __cmp__(self, other):
        assert 0, (self,other)
    def promote(self, val):
        if isinstance(val,Float):
            return val
        if isinstance(val,str):
            return VarFloat(val)
        if type(val) in (int,float,Interval):
            val = self.scalar_promote(val)
        else:
            return None
        if val==self.scalar_zero:
            return ZeroFloat()
        if val==self.scalar_one:
            return OneFloat()
        return Float(val)
    promote = classmethod(promote)
    def expr(self):
        return repr(self.val)
    def deeplen(self):
        return 1
    def uniqlen(self,cache={}):
        if self in cache:
            return 0
        else:
            cache[self]=None
            return 1
    def uniq_form(self, float_cache={}, var_cache={}, var_floats=[]):
        """
            Generate a sequence of unique "assign" statements
                float_cache: map Float -> VarFloat
                var_cache: map VarFloat -> Float  (an "assignment")
                var_floats: what order to assign the vars in
        """
        if self not in float_cache:
            var = VarFloat.temp()
            float_cache[self]=var
            var_cache[var]=self # ASSIGNMENT
            var_floats.append(var)
        return float_cache, var_cache, var_floats
#    def begin_computation(cls):
#        ctx = ComputationContext()
#        return ctx
#    begin_computation = classmethod(begin_computation)
#    def end_computation(cls):
#        pass
#        # flush cache:
#        # cls.cache = {}
#    end_computation = classmethod(end_computation)
    def get_func(self, name, free_vars=[], dumpfile=None):
        float_cache, var_cache, var_floats = self.uniq_form()
        lines = []
        lines.append( 'def %s(%s):' % (name, ', '.join( [str(v) for v in free_vars] )))
        for var in var_floats:
            lines.append( "  %s = %s" % (var, var_cache[var]) )
        lines.append( "  return %s" % (var_floats[-1]) )
        src = '\n'.join(lines)
        if dumpfile is not None:
            print >>dumpfile, src
        exec src
        return locals()[name]
    def rewrite(self, var, newvar):
        return self
    def __getattr__(self, name):
        "eg. x.sin() => SinFloat(x) "
        if name.startswith('__'):
            raise AttributeError
        orig_name = name
        name = name[0].upper()+name[1:]+'Float' # eg. SinFloat
        cls = Float.unop_cache.get(name)
        if cls is None:
            cls = type( name, (PostfixOpFloat,), {'op':orig_name} )
            Float.unop_cache[name]=cls
        return lambda : cls(self)
#    def hull(self, other):
#        other = Float.promote(other)
#        if other is None:
#            raise NotImplementedError
#        return HullFloat(self, other)
    def __str__(self):
        return self.expr()
    def __pos__( self ):
        return self
    def __neg__( self ):
        return NegFloat(self)
    def __add__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return self
        return AddFloat(self, other)
    def __radd__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return self
        return AddFloat(other, self)
    def __sub__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return self
        return SubFloat(self, other)
    def __rsub__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return -self
        return SubFloat(other, self)
    def __mul__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return ZeroFloat()
        return MulFloat(self, other)
    def __rmul__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return ZeroFloat()
        return MulFloat(other, self)
    def __div__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        if type(other) == ZeroFloat:
            return ZeroFloat()
        if type(other) == OneFloat:
            return self
        return DivFloat(self, other)
    def __rdiv__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return DivFloat(other, self)
    def __pow__( self, other ):
        if type(other) != int:
            other = Float.promote(other)
            if other is None:
                return NotImplemented
            return PowFloat(self, other)
        return IntPowFloat(self, other)

class ZeroFloat(Float):
    _zero_float = None
    def __new__(cls):
        if cls._zero_float is None:
            cls._zero_float = Float.__new__(ZeroFloat,cls.scalar_zero)
#        assert type(ZeroFloat._zero_float) == ZeroFloat, type(ZeroFloat._zero_float)
        return cls._zero_float
#    def __init__(self):
#        pass
    def __neg__( self ):
        return self
    def __add__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other
    def __radd__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other
    def __sub__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return -other
    def __rsub__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other
    def __mul__( self, other ):
        return self
    def __rmul__( self, other ):
        return self
    def __div__( self, other ):
        return self
    def __rdiv__( self, other ):
        raise ZeroDivisionError


class OneFloat(Float):
    _one_float = None
    def __new__(cls):
        if cls._one_float is None:
            cls._one_float = Float.__new__(OneFloat,cls.scalar_one)
        return cls._one_float
#    def __init__(self):
#        pass
    def __mul__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other
    def __rmul__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other
    def __rdiv__( self, other ):
        other = Float.promote(other)
        if other is None:
            return NotImplemented
        return other


class VarFloat(Float):
    uid = 0
    def __new__(cls, name):
        ob = object.__new__(cls)
        ob._hash = None
        ob.name = name
        if ob in cls.cache:
            return cls.cache[ob]
        cls.cache[ob] = ob
        return ob
#    def __init__(self, name):
##        self.name = name
#        pass
    def temp(cls):
        cls.uid += 1
        return cls('temp_%d'%cls.uid)
    temp = classmethod(temp)
    def hash(self):
        return hash(self.name)
    def __eq__(self, other):
        return type(self)==type(other) and self.name==other.name
#    def __str__(self):
#        return "VarFloat('%s')"%self.name
#    __repr__ = __str__
    def __repr__(self):
        return str(self)
    def expr(self):
        return self.name
    def rewrite(self, var, newvar):
        if var==self:
            return newvar
        return self

class UnaryOpFloat(Float):
    def __new__(cls, a):
        ob = object.__new__(cls)
        ob._hash = None
        ob.a = a
        if ob in cls.cache:
            return cls.cache[ob]
        cls.cache[ob] = ob
        return ob
#    def __init__(self, a):
#        assert isinstance(a,Float)
    def hash(self):
        return hash((self.op,self.a))
    def __eq__(self, other):
        return type(self)==type(other) and \
            (self.op,self.a)==(other.op,other.a)
    def expr(self):
        return "(%s(%s))"%(self.op, self.a.expr())
    def deeplen(self):
        return 1+self.a.deeplen()
    def uniqlen(self,cache={}):
        if self in cache:
            return 0
        else:
            cache[self]=self
            return 1+self.a.uniqlen(cache)
    def uniq_form(self, float_cache={}, var_cache={}, var_floats=[]):
        if self not in float_cache:
            self.a.uniq_form(float_cache,var_cache,var_floats)
            var = VarFloat.temp()
            a = float_cache[self.a]
            float_cache[self]=var
            var_cache[var]=self.__class__(a) # ASSIGNMENT
            var_floats.append(var)
        return float_cache, var_cache, var_floats
    def rewrite(self, var, newvar):
        return self.__class__( self.a.rewrite( var, newvar ) )

class PostfixOpFloat(UnaryOpFloat):
    def expr(self):
        return "%s.%s()"%(self.a.expr(), self.op)


class NegFloat(UnaryOpFloat):
    op = '-'
# we also build other subclasses of UnaryOpFloat
# dynamically, eg. SinFloat

class BinOpFloat(Float):
    def __new__(cls, a, b):
        ob = object.__new__(cls)
        ob._hash = None
        ob.a = a
        ob.b = b
        if ob in cls.cache:
            return cls.cache[ob]
        cls.cache[ob] = ob
        return ob
#    def __init__(self, a, b):
#        assert isinstance(a,Float), repr(a)
#        assert isinstance(b,Float), repr(b)
    def hash(self):
        return hash((self.op,self.a,self.b))
    def __eq__(self, other):
        return type(self)==type(other) and \
            (self.op,self.a,self.b)==(other.op,other.a,other.b)
    def expr(self):
        return '(%s %s %s)' % (self.a.expr(), self.op, self.b.expr())
    def deeplen(self):
        return 1+self.a.deeplen()+self.b.deeplen()
    def uniqlen(self,cache={}):
        if self in cache:
            return 0
        else:
            cache[self]=self
            return 1+self.a.uniqlen(cache)+self.b.uniqlen(cache)
    def uniq_form(self, float_cache={}, var_cache={}, var_floats=[]):
        if self not in float_cache:
            self.a.uniq_form(float_cache,var_cache,var_floats)
            self.b.uniq_form(float_cache,var_cache,var_floats)
            var = VarFloat.temp()
            a = float_cache[self.a]
            b = float_cache[self.b]
            float_cache[self]=var
            var_cache[var]=self.__class__(a,b) # ASSIGNMENT
            var_floats.append(var)
        return float_cache, var_cache, var_floats
    def rewrite(self, var, newvar):
        return self.__class__(
            self.a.rewrite( var, newvar ),
            self.b.rewrite( var, newvar ),
        )

class AddFloat(BinOpFloat):
    op = '+'
class SubFloat(BinOpFloat):
    op = '-'
class MulFloat(BinOpFloat):
    op = '*'
class DivFloat(BinOpFloat):
    op = '/'
class PowFloat(BinOpFloat):
    op = '**'
class IntPowFloat(BinOpFloat):
    op = '**'
#    def __init__(self, a, b):
#        assert isinstance(a,Float), repr(a)
#        assert isinstance(b,int), repr(b) # <--- second arg is an int <---
    def expr(self):
        return '(%s %s %s)' % (self.a.expr(), self.op, self.b)
    def deeplen(self):
        return 1+self.a.deeplen()+1
    def uniqlen(self,cache={}):
        if self in cache:
            return 0
        else:
            cache[self]=self
            return 1+self.a.uniqlen(cache)
    def uniq_form(self, float_cache={}, var_cache={}, var_floats=[]):
        if self not in float_cache:
            self.a.uniq_form(float_cache,var_cache,var_floats)
            var = VarFloat.temp()
            a = float_cache[self.a]
            b = self.b
            float_cache[self]=var
            var_cache[var]=self.__class__(a,b) # ASSIGNMENT
            var_floats.append(var)
        return float_cache, var_cache, var_floats
    def rewrite(self, var, newvar):
        return self.__class__( self.a.rewrite( var, newvar ), self.b )

#class PostfixBinOpFloat(BinOpFloat):
#    def expr(self):
#        return "%s.%s(%s)"%(self.a.expr(), self.op, self.b.expr())
#class HullFloat(PostfixBinOpFloat):
#    op = 'hull'

