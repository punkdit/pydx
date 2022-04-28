#!/usr/bin/env python

#    mjet.py : multivariate automatic differentiation module.
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

from pydx.scalar import fmath
from pydx.scalar import Scalar, set_interval_scalar, set_symbolic_scalar, restore_scalar
from pydx.scalar.symbolic import Float

try:
    import gmpy
    from gmpy import mpf
except ImportError:
    print " ** no gmpy found ** "

try:
    from scalar.xinterval import Interval, get_default_prec, set_default_prec, DisorderException
except ImportError, e:
    print " ** could not import xinterval ** "
    print e

__module__ = sys.modules[__name__]

Scalar.clients.append(__module__)

#
###############################################################################
#

#class slice(slice):
#    def __hash__(self):
#        return hash((self.start,self.stop,self.step))
# TypeError: Error when calling the metaclass bases
#         type 'slice' is not an acceptable base type

# do these two functions belong in this module ??
def xcross(idxs):
    if len(idxs):
        idx = idxs[0]
        if type(idx)==slice:
            assert idx.step is None, "Not Implemented"
            for idx in range(idx.start, idx.stop):
                for rest in xcross(idxs[1:]):
                    yield (idx,)+rest
        else:
            for rest in xcross(idxs[1:]):
                yield (idx,)+rest
    else:
        yield ()

def cross(idxs):
    if len(idxs):
        for idx in range(0,idxs[0]):
            for rest in cross(idxs[1:]):
                yield (idx,)+rest
    else:
        yield ()

### This cache doesn't seem to speed things up any.
#_multi_range_order_cache = {}
#def multi_range_order(*args):
#    idxss = _multi_range_order_cache.get(args)
#    if idxss is None:
#        idxss = tuple(idxs for idxs in _multi_range_order(*args))
#        _multi_range_order_cache[args]=idxss
#    return idxss
#def _multi_range_order(rank,start_order,stop_order):
def multi_range_order(rank,start_order,stop_order):
    """
        yield all multi indexes with sum from start_order to stop_order inclusive.
    """
    assert 0<=start_order<=stop_order
    if rank==0:
        assert start_order==stop_order==0
        yield ()
    elif rank==1:
        for i in range(start_order,stop_order+1):
            yield (i,)
    else:
        for i in range(0,stop_order+1):
            for rest in multi_range_order(rank-1,max(0,start_order-i),stop_order-i):
                yield (i,)+rest

def multi_range(idxs):
    if len(idxs):
        for idx in range(0,idxs[0]+1):
            for rest in multi_range(idxs[1:]):
                yield (idx,)+rest
    else:
        yield ()

def multi_zero(rank):
    return (0,)*rank
def multi_unit(rank,i):
    idxs = [0]*rank
    idxs[i] = 1
    return tuple(idxs)
def multi_unit_under(idxs):
    " return unit _idxs s.t. _idxs<=idxs "
    i = 0
    while idxs[i]==0: # IndexError means that idxs all 0
        i += 1
    _idxs=[0]*len(idxs)
    _idxs[i]=1
    return tuple(_idxs)

def multi_lexcmp(a_idxs, b_idxs):
    " lexicographical ordering "
    assert len(a_idxs)==len(b_idxs)
    i = 0
    while i < len(a_idxs) and a_idxs[i] == b_idxs[i]:
        i += 1
    if i==len(a_idxs):
        return 0
    return cmp(a_idxs[i],b_idxs[i])

def multi_cmp(a_idxs, b_idxs):
    le = True
    ge = True
    eq = True
    for a_idx,b_idx in zip(a_idxs,b_idxs):
        if a_idx>b_idx:
            eq = False
            le = False
        elif a_idx<b_idx:
            eq = False
            ge = False
    if eq:
        return 0
    if le:
        return -1
    if ge:
        return 1

def multi_add(a_idxs, b_idxs):
    assert len(a_idxs)==len(b_idxs)
    return tuple( a_idxs[i]+b_idxs[i] for i in range(len(a_idxs)) )

def multi_sub(a_idxs, b_idxs):
    assert len(a_idxs)==len(b_idxs)
    c_idxs = tuple( a_idxs[i]-b_idxs[i] for i in range(len(a_idxs)) )
    assert multi_cmp(c_idxs,multi_zero(len(a_idxs))) >= 0
    return c_idxs

def multi_partitions(idxs):
    rank = len(idxs)
    zero = multi_zero(rank)
    if idxs == zero:
        yield []
    else:
        for _idxs in multi_range( idxs ):
            if _idxs == zero:
                continue
            for part in multi_partitions( multi_sub( idxs, _idxs ) ):
                if part and multi_lexcmp( _idxs, part[0] ) > 0: # XX use extra parameter
                    continue
                yield [_idxs]+part

def n_choice(n,k):
    " n_choice(n,k): n choose k: n!/( k!(n-k)! ) "
    assert k<=n
    if k==0:
        return 1
    if k==1:
        return n
    i = 1
    for j in range(n-k+1,n+1):
        i *= j
    for j in range(1,k+1):
        i /= j
    return i


class MJet(object):
    scalar_type = float
    scalar_promote = float
    scalar_zero = 0.0
    scalar_one = 1.0
    def __init__( self, rank ):
        self.rank = rank
    def promote( self, item ):
        if not isinstance(item,MJet):
#            assert type(item)==self.scalar_type, (repr(item),type(item))
#            assert isinstance(item,self.scalar_type), (repr(item),type(item))
            item = self.scalar_promote(item)
            x = Jet({(0,)*self.rank:item})
            return x
        elif item.rank == self.rank:
            return item
        assert self.rank > item.rank
#        assert isinstance(item,Jet), "lazy promote not implemented (item.rank=%d)"%item.rank
        if isinstance(item,Jet):
            x = Jet(rank = self.rank)
            extra_zeros = (0,)*(self.rank-item.rank)
            for idxs, value in item.cs.items():
                x[ idxs+extra_zeros ] = value
            return x
        else:
            return PromotedJet( item, self.rank )
    def str(self, n):
        return "<%s>" % ' '.join(str(self[i]) for i in range(n) )
    def repr(self, n):
        return "<%s>" % ' '.join(repr(self[i]) for i in range(n) )
    def get_scalar_component(self):
        return self[(0,)*self.rank]
    scalar_component=property(get_scalar_component)

    def bdot(cls, a, b, n, i):
        """
            @param a: MJet
            @param b: MJet
            @param n: tuple
            @param i: int
        """
        rank = len(n)
        delta_i = multi_unit(rank, i)
        assert n[i]>0
        r = sum(( (n[i]-j[i])*a[j]*b[multi_sub(n,j)] for j in multi_range(multi_sub(n,delta_i))),
                        cls.scalar_zero )
        r = r / n[i]
        return r
    bdot = classmethod(bdot)

#    def bdot(cls, P, m, Q, k ):
#        """
#            @param P: MJet
#            @param m: tuple
#            @param Q: MJet
#            @param k: tuple
#        """
#        rank = len(m)
#        assert multi_cmp(m,k)<=0
#        
#        return r
#    bdot = classmethod(bdot)

    def expand( self, x, order ):
        assert len(x)==self.rank
        r = self.scalar_zero
        for j in multi_range_order(self.rank, 0, order):
            s = self.scalar_one
            for i in range(self.rank):
                if j[i]>0:
                    s *= x[i]**j[i]
            r = r + self[j] * s # XX FACTORIZE THIS XX
        return r

    #
    # On order=5 this has around 0.001% more error: (???)
    def expand_factorized( self, x, order ):
        # factored taylor series
        assert len(x)==self.rank
        zero_idxs = (0,)*self.rank
        
        terms = {}
        for j in multi_range_order(self.rank, order, order):
            terms[j] = self[j]
        while order:
            order -= 1
            _terms = {} # new terms for next iteration
            for j in multi_range_order(self.rank, order, order):
                _terms[j] = self[j]
            for idxs in terms:
                i = 0
                while idxs[i]==0:
                    i += 1
                value = terms[idxs] * x[i]
                _idxs = multi_sub(idxs,multi_unit(self.rank,i))
                _terms[_idxs] = _terms[_idxs] + value
            terms = _terms
        assert len(terms)==1
        return terms[zero_idxs]

    def expand_err( self, x, order, err ):
        assert len(x)==self.rank
#        assert len(err)==order+2 # what should this be ?
        assert self.scalar_type in ( Interval, Float )
#        print "expand_err self:", ' '.join([str(self[i]) for i in range(3)])
#        print "expand_err x=%s order=%d err=%s" % (x[0], order, err[0])
        r = self.scalar_zero
        for j in multi_range_order(self.rank, 0, order):
            s = self.scalar_one
            for i in range(self.rank):
                if j[i]>0:
                    s *= x[i]**j[i]
            r = r + self[j] * s
#            print "expand_err j=%s r=%s" % (j,r)
        for i,j in enumerate(multi_range_order(self.rank, order+1, order+1)):
            s = self.scalar_one
            for k in range(self.rank):
                if j[k]>0:
                    s *= x[k]**j[k]
            r = r + err[i] * s
#            print "expand_err r=%s" % r
        assert i==len(err)-1
        return r

    def expand_err_factorized( self, x, order, err ):
        # factored taylor series (less accurate)
        assert len(x)==self.rank
#        assert len(err)==order+2 # what should this be ?
        zero_idxs = (0,)*self.rank
        
        order += 1
        terms = {}
        for i,j in enumerate(multi_range_order(self.rank, order, order)):
            terms[j] = err[i]
        while order:
            order -= 1
            _terms = {} # new terms for next iteration
            for j in multi_range_order(self.rank, order, order):
                _terms[j] = self[j]
            for idxs in terms:
                i = 0
                while idxs[i]==0:
                    i += 1
                value = terms[idxs] * x[i]
                _idxs = multi_sub(idxs,multi_unit(self.rank,i))
                _terms[_idxs] = _terms[_idxs] + value
            terms = _terms
        assert len(terms)==1
        return terms[zero_idxs]

    def __eq__(self, other):
        return repr(self)==repr(other)
    def __pos__( self ):
        return self
    def __neg__( self ):
        return NegJet(self)
    def __add__( self, other ):
        return AddJet(self, other)
    def __radd__( self, other ):
        other = self.promote(other)
        if other is None:
            return NotImplemented
        return AddJet(other, self)
    def __sub__( self, other ):
        return SubJet(self, other)
    def __rsub__( self, other ):
        other = self.promote(other)
        if other is None:
            return NotImplemented
        return SubJet(other, self)
    def __mul__( self, other ):
        return MulJet(self, other)
    def __rmul__( self, other ):
        other = self.promote(other)
        if other is None:
            return NotImplemented
        return MulJet(other, self)
    def __div__( self, other ):
        return DivJet(self, other)
    def __rdiv__( self, other ):
        other = self.promote(other)
        if other is None:
            return NotImplemented
        return DivJet(other, self)
    def __pow__( self, other ):
        if other==0:
            r = Jet(rank=self.rank)
            r[(0,)*self.rank] = self.scalar_one
        elif other == 1:
            r = self
        elif other == 2:
            r = SqrJet(self)
        elif int(other)==other:
            return IntPowJet(self, other)
        elif other==0.5:
            r = SqrtJet(self)
        else:
            other = self.scalar_promote(other)
            r = (other*self.log()).exp()
        return r
    def exp(self):
        return ExpJet(self)
    def log(self):
        return LogJet(self)
    def sqrt(self):
        return SqrtJet(self)
    def atan(self):
        return ATanJet(self)
    def asin(self):
        return ASinJet(self)
    def acos(self):
        return ACosJet(self)
    def sin(self):
        return SinJet(self)
    def cos(self):
        return CosJet(self)
    def tan(self):
        return SinJet(self)/CosJet(self)
    def cot(self):
        return CosJet(self)/SinJet(self)
    def reciprocal(self):
        return self.scalar_one/self
    recip = reciprocal
    def sec(self):
        return self.cos().reciprocal()
    def sinh(self):
        return (self.exp()-(-self).exp())/(self.scalar_one*2)
    def cosh(self):
        return (self.exp()+(-self).exp())/(self.scalar_one*2)
    def tanh(self):
        return (self.exp()-(-self).exp()) / (self.exp()+(-self).exp())
Scalar.clients.append( MJet )

class JetOp(MJet): # Lazy Jet with item cache
    def __init__( self, rank ):
        MJet.__init__( self, rank )
        self._cache = {} # speeds up things by *30
    def __getitem__( self, idxs ):
        if type(idxs)==int:
            idxs=idxs,
        assert type(idxs)==tuple
        assert len(idxs)==self.rank
        item = self._cache.get(str(idxs)) # XX slice objects are not hashable :(
        if item is None:
            if [idx for idx in idxs if type(idx)==slice]:
                item = SliceJet(self, idxs)
            else:
                item = self.getitem(idxs)
            self._cache[str(idxs)] = item
        return item
    def __str__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,self.a,self.b)
    def __repr__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,repr(self.a),repr(self.b))
    __setitem__ = None

class PromotedJet(JetOp):
    def __init__( self, a, rank ):
        JetOp.__init__(self, rank)
        self.a = a
        assert rank > self.a.rank
    def getitem( self, idxs ):
        if not idxs[self.a.rank:self.rank]==(0,)*(self.rank-self.a.rank):
            # and it hasn't been set, so:
            return self.scalar_zero
        return self.a[ idxs[:self.a.rank] ]
    def __setitem__( self, idxs, value ):
        assert type(idxs)==tuple
        self._cache[str(idxs)] = value

class AddJet(JetOp):
    def __init__( self, a, b ):
        JetOp.__init__(self, a.rank)
        b = self.promote(b)
        assert b.rank == a.rank
        self.a = a
        self.b = b
    def getitem( self, idxs ):
        assert type(idxs)==tuple
        assert isinstance(self.a[idxs],self.scalar_type), (self.a[idxs],self.scalar_type)
        assert isinstance(self.b[idxs],self.scalar_type)
        return self.a[idxs]+self.b[idxs]

class SubJet(JetOp):
    def __init__( self, a, b ):
        JetOp.__init__(self, a.rank)
        b = self.promote(b)
        assert b.rank == a.rank
        self.a = a
        self.b = b
    def getitem( self, idxs ):
        assert type(idxs)==tuple
        assert isinstance(self.a[idxs],self.scalar_type)
        assert isinstance(self.b[idxs],self.scalar_type)
        return self.a[idxs]-self.b[idxs]

class NegJet(JetOp):
    def __init__( self, a ):
        JetOp.__init__(self, a.rank)
        self.a = a
    def getitem( self, idxs ):
        assert type(idxs)==tuple
        assert isinstance(self.a[idxs],self.scalar_type)
        return -self.a[idxs]
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.a)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.a))

class MulJet(JetOp):
    def __init__( self, a, b ):
        JetOp.__init__(self, a.rank)
        b = self.promote(b)
        assert b.rank == a.rank
        self.a = a
        self.b = b
    def getitem( self, idxs ):
        assert type(idxs)==tuple
        idxss = [ _ for _ in multi_range(idxs) ]
#        _idxss = [ tuple([ idxs[i]-_idxs[i] for i in range(len(idxs)) ]) for _idxs in idxss ]
        _idxss = [ multi_sub(idxs,_idxs) for _idxs in idxss ]
        return sum( [ self.a[a_idxs] * self.b[b_idxs] for a_idxs, b_idxs in zip(_idxss,idxss) ], self.scalar_zero )

class DivJet(JetOp):
    def __init__( self, b, c ):
        JetOp.__init__(self, b.rank)
        c = self.promote(c)
        assert c.rank == b.rank
        self.b = b
        self.c = c
    def __str__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,self.b,self.c)
    def __repr__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,repr(self.b),repr(self.c))
    def getitem( self, idxs ):
        assert type(idxs)==tuple
        idxss = [ _ for _ in multi_range(idxs) if _ != (0,)*self.rank ]
#        _idxss = [ tuple([ idxs[i]-_idxs[i] for i in range(len(idxs)) ]) for _idxs in idxss ]
        _idxss = [ multi_sub(idxs,_idxs) for _idxs in idxss ]
        r = self.b[idxs]
        r = r - sum( [ self.c[c_idxs] * self[self_idxs] 
            for c_idxs, self_idxs in zip(idxss,_idxss) ], self.scalar_zero )
        r = r / self.c[ (0,)*self.rank ]
        return r

class IntPowJet(JetOp):
    def __init__( self, b, n ):
        assert type(n)==int
        JetOp.__init__(self, b.rank)
        self.b = b
        self.n = n
    def __str__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,self.b,self.n)
    def __repr__(self):
        return "%s(%s,%s)"%(self.__class__.__name__,repr(self.b),repr(self.n))
    def getitem( self, idxs ):
        parts = multi_partitions(idxs)
        r = self.scalar_zero
        mzero = multi_zero(self.rank)
        for part in parts:
            if len(part)>self.n:
                # only use up to n factors. XX add parameter to multi_partitions XX
                continue
            _r = self.scalar_one
            if len(part)<self.n:
                # make up the extra factors
                _r = self.b[ mzero ]**(self.n-len(part))
            last_idxs = None
            exponents = {}
            for idxs in part:
                exponents[idxs] = exponents.get(idxs,0) + 1
            n = self.n # how many to choose from
            for idxs, exponent in exponents.items():
                _r = n_choice(n, exponent) * _r * ( self.b[idxs] ** exponent )
                n -= exponent # we chose this many
            r = r + _r
        return r

#class PowJet(JetOp):
#    def __init__( self, other, alpha ):
#        JetOp.__init__(self, other.rank)
#        self.other = other
#        self.alpha = self.scalar_promote(alpha)
#    def __str__(self):
#        return "%s(%s,%s)"%(self.__class__.__name__,self.other,self.alpha)
#    def __repr__(self):
#        return "%s(%s,%s)"%(self.__class__.__name__,repr(self.b),repr(self.n))
#    def getitem( self, n ):
#        raise NotImplementedError
#        zero_idxs = multi_zero(self.rank)
#        alpha = self.alpha
#        a, b = self, self.other
#        if n==zero_idxs:
#            return alpha * (b[zero_idxs]**(alpha-self.scalar_one))
#
#        return r

class UnaryJetOp(JetOp):
    def __init__( self, b ):
        JetOp.__init__(self, b.rank)
        self.b = b
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.b)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.b))

class ExpJet(UnaryJetOp):
    """
        a = exp(b)
    """
    def getitem( self, n ):
        zero_idxs = multi_zero(self.rank)
        a, b = self, self.b
        if n==zero_idxs:
            return fmath.exp(b[n])
        i = 0
        while n[i]==0:
            i += 1
        r = MJet.bdot( a, b, n, i )
        return r

class LogJet(UnaryJetOp):
    def getitem( self, k ):
        # h = log(u)
        zero_idxs = multi_zero(self.rank)
        h, u = self, self.b
        if k==zero_idxs:
            return fmath.log(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            if j != zero_idxs and k[i]>j[i]:
                r = r + (k[i]-j[i])*u[j]*h[multi_sub(k,j)]
        r = r/k[i]
        r = (u[k] - r) / u[zero_idxs]
        return r
        
class SqrtJet(UnaryJetOp):
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        h, u = self, self.b
        if k==zero_idxs:
            return fmath.sqrt(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            if j != zero_idxs and k[i]>j[i]:
                r = r + (k[i]-j[i])*h[j]*h[multi_sub(k,j)]
        r = r/k[i]
        r = (u[k]/(2*self.scalar_one) - r) / h[zero_idxs]
        return r
        
class SqrJet(UnaryJetOp):
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        a, b = self, self.b
        if k==zero_idxs:
            return fmath.sqr(b[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        ei = multi_unit(self.rank,i)
        ki = multi_sub(k,ei)
        for j in multi_range( ki ):
            r = r + (j[i]+1.0) * b[ multi_sub(ki,j) ] * b[ multi_add(ei,j) ]
        return 2.0 * r / k[i]


class ATanJet(UnaryJetOp):
    def __init__( self, u ):
        JetOp.__init__(self, u.rank)
        self.u = u
        self.v = self.scalar_one / ( self.scalar_one + u*u )
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.u)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.u))
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        h, u, v = self, self.u, self.v
        if k==zero_idxs:
            return fmath.atan(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            r = r + (k[i]-j[i])*v[j]*u[multi_sub(k,j)]
        r = r/k[i]
        return r

class ASinJet(UnaryJetOp):
    def __init__( self, u ):
        JetOp.__init__(self, u.rank)
        self.u = u
        self.v = self.scalar_one / ( self.scalar_one - u*u ).sqrt()
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.u)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.u))
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        h, u, v = self, self.u, self.v
        if k==zero_idxs:
            return fmath.asin(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            r = r + (k[i]-j[i])*v[j]*u[multi_sub(k,j)]
        r = r/k[i]
        return r

class ACosJet(UnaryJetOp):
    def __init__( self, u ):
        JetOp.__init__(self, u.rank)
        self.u = u
        raise NotImplementedError
        self.v = self.scalar_one / ( self.scalar_one - u*u ).sqrt()
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.u)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.u))
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        h, u, v = self, self.u, self.v
        if k==zero_idxs:
            return fmath.asin(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            r = r + (k[i]-j[i])*v[j]*u[multi_sub(k,j)]
        r = r/k[i]
        return r

class SinJet(UnaryJetOp):
    def __init__( self, u, cosu=None ):
        JetOp.__init__(self, u.rank)
        self.u = u
        if cosu is None:
            cosu = CosJet(u,self)
        self.cosu = cosu
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.u)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.u))
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        sinu, u, cosu = self, self.u, self.cosu
        if k==zero_idxs:
            return fmath.sin(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            r = r + (k[i]-j[i])*cosu[j]*u[multi_sub(k,j)]
        r = r/k[i]
        return r

class CosJet(UnaryJetOp):
    def __init__( self, u, sinu=None ):
        JetOp.__init__(self, u.rank)
        self.u = u
        if sinu is None:
            sinu = SinJet(u,self)
        self.sinu = sinu
    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__,self.u)
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,repr(self.u))
    def getitem( self, k ):
        zero_idxs = multi_zero(self.rank)
        cosu, u, sinu = self, self.u, self.sinu
        if k==zero_idxs:
            return fmath.cos(u[k])
        i = 0
        while k[i]==0:
            i += 1
        r = self.scalar_zero
        for j in multi_range( multi_sub(k,multi_unit(self.rank,i)) ):
            r = r - (k[i]-j[i])*sinu[j]*u[multi_sub(k,j)]
        r = r/k[i]
        return r

class SliceJet(JetOp):
    """
        SliceJet(a,idxs): take a (lazy) slice out of 'a'
    """
    def __init__( self, a, idxs ):
        assert type(idxs)==tuple
        # every slice gives us a rank
        rank = len([idx for idx in idxs if type(idx)==slice])
        JetOp.__init__( self, rank=rank )
        self.a = a
        self.full_idxs = idxs
        assert len(idxs)==a.rank # full index into a
        idx_map = {}
        i = 0
        for j, idx in enumerate(idxs):
            if type(idx)==slice:
                idx_map[i] = j
                i += 1
                assert idx==slice(None), " shift not implemented "
        self.idx_map = idx_map
    def getitem( self, idxs ):
        full_idxs = list(self.full_idxs) # copy
        for i, idx in enumerate(idxs):
            full_idxs[self.idx_map[i]] = idx
        return self.a[tuple(full_idxs)]
    def __str__(self):
        return "%s(%s, rank=%d, idxs=%s)"%(self.__class__.__name__,self.a,self.rank,self.full_idxs)
    def __repr__(self):
        return "%s(%s, rank=%d, idxs=%s)"%(self.__class__.__name__,repr(self.a),self.rank,self.full_idxs)
        
class SprayItem(JetOp):
    " Created by Spray. These are the (MJet) components of a Spray object. "
    def __init__(self, spray, idx, yi):
        JetOp.__init__(self, rank=1) # yi.rank+1 ?
        self.spray = spray
        self.idx = idx
        self.yi = MJet(1).promote(yi)
    def getitem(self, order):
#        assert type(order)==int, "not implemented"
        #print "SprayItem.getitem: order =", order
        order, = order
        if order==0:
            return self.yi[(0,)]
        dyi = self.spray.dy[self.idx]
        dyi = MJet(1).promote(dyi)
        r = (self.scalar_one/order) * dyi[order-1]
        #print "SprayItem.getitem r =", r
        return r
class Spray(object):
    " Behaves like a list (vector) of MJet's "
    " Defined by an ODE: dy=f(x,y)"
    def __init__(self, f, x, y):
        self.f = f
        self.x = x
        self.y = [SprayItem(self,idx,yi) for idx,yi in enumerate(y)]
        self.dy = f( x, self.y ) # a list of MJet's
    def __getitem__(self, idx):
        return self.y[idx]

class Jet(MJet):
    """
        This is a concrete Jet. All components are stored in a dictionary.
    """
    def __init__( self, cs={}, rank=None ):
        if type(cs) in (list,tuple):
            if cs and type(cs[0]) in (list,tuple):
                raise NotImplementedError
            cs = dict( ((i,),c) for i,c in enumerate(cs) )
        if type(cs)!=dict:
            raise ValueError, "expected dict, list or tuple, but got %s"%cs
        if rank is None and not cs:
            raise ValueError, "can't determine rank"
        elif rank is None:
            rank = len(cs.keys()[0])
        for key in cs.keys():
            assert len(key)==rank
        MJet.__init__(self, rank)
        self.bound = [0]*rank
        self.cs = {} # coefficients, components, ...
        for idxs, value in cs.items():
            # assert type(value) == self.scalar_type ??
            self[idxs] = self.scalar_promote(value)
        
    def __str__(self):
        ss = []
        for idxs in xcross( [slice(0,n+1) for n in self.bound] ):
            ss.append( '%s:%s' % (''.join([str(idx) for idx in idxs]), self[idxs] ) )
        return '<%s>'%' '.join(ss)
    def __repr__(self):
        ss = []
        for idxs in xcross( [slice(0,n+1) for n in self.bound] ):
            ss.append( '%s:%s' % (''.join([str(idx) for idx in idxs]), repr(self[idxs]) ) )
        return '<%s>'%' '.join(ss)
    def __getitem__( self, idxs ):
        if type(idxs)==int:
            idxs=idxs,
        assert type(idxs)==tuple
        assert len(idxs)==self.rank
        # every slice gives us a rank
        rank = len([idx for idx in idxs if type(idx)==slice])
        if rank:
            # Non-scalar return type
            idx_map = {}
            idxs = list(idxs)
            slice_i = []
            for i, idx in enumerate(idxs):
                if type(idx)==slice:
                    slice_i.append(i)
                    assert idx==slice(None), "shift not implemented"
                    idx = slice(0,self.bound[i]+1,None)
                    idxs[i] = idx
            x = Jet( rank = rank )
            for full_idxs in xcross(idxs):
                idx = [full_idxs[i] for i in slice_i]
                x[tuple(idx)] = self[full_idxs]
            return x
        return self.cs.get( idxs, self.scalar_zero )
    def __setitem__( self, idxs, c ):
        if type(idxs)==int:
            idxs=idxs,
        assert type(idxs)==tuple
        assert len(idxs)==self.rank
#        assert type(c)==float
        assert isinstance(c,self.scalar_type), repr(c)
        self.cs[idxs] = c
        for i, idx in enumerate(idxs):
            self.bound[i] = max(idx, self.bound[i])
    def __cmp__(self, other):
        assert type(other)==Jet, "Not implemented"
        for idxs in self.cs.keys()+other.cs.keys():
            r = cmp( self.cs.get(idxs,0.0), other.cs.get(idxs,0.0) )
            if r != 0:
                return -1
        return 0

class Derivation(object):
    """
        given f:R^n->R
        Derivation(f):R^n->R^n
            a list of the partial derivative's of f
            extends to arbitrary Jet arguments
    """
    def __init__( self, f ):
        self.f = f
    def __call__( self, *y ):
        f = self.f
        dim = len(y)
        assert dim
        if not isinstance(y[0],MJet):
            # must be scalar .. ie. a rank-0 Jet
            y = [ MJet(0).promote(_y) for _y in y ]
        # now we need dim extra ranks
        rank0 = y[0].rank
        rank = rank0 + dim
        yy = []
        for i,_y in enumerate(y):
            _y = MJet(rank).promote(_y) # XX not lazy yet XX
            idxs = [0]*rank
            idxs[rank0 + i] = 1
            _y[ tuple(idxs) ] = 1.0 # self.scalar_one XX
            yy.append(_y)
    
        ##############
        g = f( *yy ) # 
        ##############
    
    #  print g[0,0], g[0,1], g[1,0], g[1,1]
        g = MJet(rank).promote(g)
        # now we drop back to rank0
        gs = []
        for i in range(dim):
            idxs = [0]*dim
            idxs[i] = 1
            idxs = tuple([slice(None)]*rank0+idxs)
            g0 = g[idxs]
            g0 = MJet(rank0).promote(g0)
            gs.append(g0)
        return gs
            
#
###############################################################################
#
