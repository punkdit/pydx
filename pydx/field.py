#!/usr/bin/env python

#    field.py: tensor field objects.
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

from pydx import scalar
from pydx.tensor import Tensor, up, dn
from pydx.mjet import MJet, cross 
from pydx.scalar import set_symbolic_scalar, restore_scalar
from pydx.scalar.symbolic import VarFloat, ComputationContext

class TensorField(object):
    """
        Instances are callable, returning a Tensor object.
        Abstract base class.
    """
    scalar_zero = 0.0
    scalar_one = 1.0
    scalar_promote = float
    scalar_type = float
    up = 1; dn = -1;
    def __init__( self, valence = (), dim = 4, g=None ):
        """
            @param g: metric
        """
        # for now, all Tensors are "square"
        self.valence = tuple(valence) # sequence of up/dn
        for updn in self.valence:
            assert updn in (self.up,self.dn), updn
        self.rank = len(valence)
        self.shape = (dim,)*self.rank
        self.dim = dim
        assert g is None or isinstance(g, TensorField)
        self.g = g
        if g is not None and 'uu' not in g.__dict__:
            g.uu = None 
        valence_attr = ''.join([{Tensor.up:'u',Tensor.dn:'d'}[v] for v in valence])
        self.__dict__[valence_attr] = self
    def __str__( self ):
        return "%s( %s, %s )"%(self.__class__.__name__, self.valence, self.dim)
    __repr__ = __str__
    def identity(cls, valence, dim ):
        tfield = ConcreteTensorField( valence, dim )
        one_func = lambda *xs: cls.scalar_one
        zero_func = lambda *xs: cls.scalar_zero
        for idxs in tfield.genidx():
            on_diag = True
            if idxs:
                i = idxs[0]
            for j in idxs:
                if j!=i:
                    on_diag = False
            if on_diag:
                tfield[idxs] = one_func
            else:
                tfield[idxs] = zero_func
        return tfield
    identity = classmethod(identity)
    def zero(cls, valence, dim ):
        tfield = ConcreteTensorField( valence, dim )
        zero_func = lambda *xs: cls.scalar_zero
        for idxs in tfield.genidx():
            tfield[idxs] = zero_func
        return tfield
    zero = classmethod(zero)
    def genidx( self ):
        return cross( (self.dim,)*self.rank )
    def __getattr__( self, valence_s ):
        " For example: riemann.dddd is the fully contravariant riemann "
        # i am too tricky for my socks
        n = valence_s.count('u')+valence_s.count('d')
        if not n==len(valence_s)==len(self.valence):
            raise AttributeError, valence_s
        assert self.g is not None
        assert self.g.uu is not None
        t = self
        for v,idx in enumerate(valence_s):
            if v=='u' and t.valence[idx]==Tensor.dn:
                t = t.mul( self.g.uu, (idx,0) )
            elif v=='d' and t.valence[idx]==Tensor.up:
                t = t.mul( self.g, (idx,0) )
        return t
    def view( self, trans ):
        return ViewTensorField( self, trans )
    def transform( self, trans ):
        return TransformTensorField( self, trans )
    def comma( self ):
        return CommaTensorField(self)
    def transpose( self, *perm ):
        return TransposeTensorField(self, perm)
    def __pos__( self ):
        return self
    def __neg__( self ):
        return NegTensorField(self)
    def __add__( self, other ):
        return AddTensorField(self, other)
#    def __radd__( self, other ):
#        other = self.promote(other)
#        return AddTensorField(other, self)
    def __sub__( self, other ):
        return SubTensorField(self, other)
#    def __rsub__( self, other ):
#        other = self.promote(other)
#        return Sub(other, self)
#    def __mul__( self, other ):
#        return ScalarMulTensorField(self, other)
    def __rmul__( self, other ):
        return ScalarMulTensorField(other, self)
    def mul( self, other, *pairs ):
        return MulTensorField( self, other, pairs )
    def contract( self, *pairs ):
        return ContractTensorField( self, pairs )
    def compile_oldversion(self):
        "the older version of the compile method. slower, but slightly more accurate (why?!?) "
        Float.scalar_zero = Interval(0.0)
        Float.scalar_one = Interval(1.0)
        Float.scalar_promote = Interval
        Float.scalar_type = Interval
        compiledump = open('compiledump.py','w')
        set_symbolic_scalar()
        tfield = ConcreteTensorField( self.valence, self.dim )
        vs = [ VarFloat( 'x%d'%i ) for i in range(self.dim) ]
        tensor = self(*vs)
        for idxs in self.genidx():
            r = tensor[idxs]
            #if not isinstance(r,MJet):
            #  r = MJet(0).promote(r)
            r = r[()]
            #if not isinstance(r,Float):
            #  r = Float.promote(r)
            name = 'func' + ''.join(str(i) for i in idxs)
            func = r.get_func( name, vs, compiledump )
            tfield[idxs] = func
        restore_scalar()
        Float.scalar_zero = 0.0
        Float.scalar_one = 1.0
        Float.scalar_promote = float
        Float.scalar_type = float
        return tfield
    def compile(self):
        compiledump = open('compiledump.py','w')
        set_symbolic_scalar()
        vs = [ VarFloat( 'x%d'%i ) for i in range(self.dim) ]
        name = 'func'
        args = vs + [ 'tensor' ]
        ctx = ComputationContext(name, args, compiledump)
        tensor = self(*vs)
        for idxs in self.genidx():
            r = tensor[idxs]
            #if not isinstance(r,MJet):
            #  r = MJet(0).promote(r)
            r = r[()]
            #if not isinstance(r,Float):
            #  r = Float.promote(r)
            ctx.assign( 'tensor[%s]'%repr(idxs), r )
        func = ctx.finalize()
        restore_scalar()
        tfield = AdaptedTensorField( self.valence, self.dim, self.g, func )
        return tfield
scalar.Scalar.clients.append(TensorField)

class AdaptedTensorField(TensorField):
    """ Used with ComputationContext. Creates tensors via ONE function call.
    """
    def __init__( self, valence = (), dim = 4, g=None, func=None ):
        TensorField.__init__( self, valence, dim, g=g )
        self.func = func
    def __call__( self, *xs ):
        tensor = Tensor( self.valence, self.dim )
        args = list(xs)+[tensor]
        self.func(*args)
        return tensor

class ConcreteTensorField(TensorField):
    """ Components are actual (scalar valued) functions
    """
    def __init__( self, valence = (), dim = 4, g=None, elems={} ):
        TensorField.__init__( self, valence, dim, g=g )
        self.elems = dict(elems) # map tuple to value
        #for idxs in cross( (self.dim,)*self.rank ):
        #  self.elems[idxs] = lambda *xs: return scalar_zero
        self.symetry = {}
    def addsymetry(self, src_idxs, tgt_idxs):
        "when src_idxs are requesed, we look-up tgt_idxs"
        assert tgt_idxs not in self.symetry, " chain a symetry ?"
        while tgt_idxs in self.symetry:
            tgt_idxs = self.symetry[tgt_idxs]
        self.symetry[src_idxs] = tgt_idxs
    def __setitem__( self, idxs, func ):
        assert callable(func)
        if type(idxs)==int:
            idxs = idxs,
        assert idxs not in self.symetry, "symetric component!"
        self.elems[idxs] = func
    def __getitem__( self, idxs ):
        if type(idxs)==int:
            idxs = idxs,
        assert idxs not in self.symetry, "symetric component!"
        func = self.elems[idxs]
        return func
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        elems = {}
        for idxs in self.genidx():
            if idxs in self.symetry: # is this any use ??
                tgt_idxs = self.symetry[idxs]
                assert tgt_idxs not in self.symetry
                if tgt_idxs not in elems:
                    elems[tgt_idxs] = self[tgt_idxs]( *xs )
                elems[idxs] = elems[tgt_idxs]
            else:
                elems[idxs] = self[idxs]( *xs )
        tensor = Tensor( self.valence, self.dim, elems=elems )
        return tensor

class NegTensorField(TensorField):
    def __init__( self, a, g=None ):
        TensorField.__init__( self, a.valence, a.dim, g=g )
        self.a = a
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        return -self.a(*xs)

class AddTensorField(TensorField):
    def __init__( self, a, b, g=None ):
        assert a.valence == b.valence
        assert a.dim == b.dim
        TensorField.__init__( self, a.valence, a.dim, g=g )
        self.a = a
        self.b = b
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        return self.a(*xs) + self.b(*xs)

class SubTensorField(TensorField):
    def __init__( self, a, b, g=None ):
        assert a.valence == b.valence
        assert a.dim == b.dim
        TensorField.__init__( self, a.valence, a.dim, g=g )
        self.a = a
        self.b = b
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        return self.a(*xs) - self.b(*xs)

class ScalarMulTensorField(TensorField):
    def __init__( self, a, b, g=None ):
        TensorField.__init__( self, b.valence, b.dim, g=g )
        self.a = a
        self.b = b
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        return self.a * self.b(*xs)

class MulTensorField(TensorField):
    def __init__( self, a, b, pairs, g=None ):
        assert a.dim == b.dim
        for i, j in pairs:
            assert a.valence[i] != b.valence[j] # pair up with dn, and vice-versa
        ai = [ aii for aii,_ in pairs ]
        bi = [ bii for _,bii in pairs ]
#        print "MulTensorField", a.valence, b.valence, pairs
        valence = []
#        print 'valence:', valence
        for i,v in enumerate(a.valence):
            if i not in ai:
                valence.append(v)
#                print 'valence:', valence
        for i,v in enumerate(b.valence):
            if i not in bi:
                valence.append(v)
#                print 'valence:', valence
#        print 'valence:', valence
        TensorField.__init__( self, tuple(valence), a.dim, g=g )
        self.a = a
        self.b = b
        self.pairs = pairs
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        t = self.a(*xs).mul( self.b(*xs), *self.pairs )
        assert t.valence == self.valence
        return t

class ContractTensorField(TensorField):
    def __init__( self, a, pairs, g=None ):
        ii = [ aii for aii,_ in pairs ] + [ bii for _,bii in pairs ]
        valence = []
        for i,v in enumerate(a.valence):
            if i not in ii:
                valence.append(v)
        TensorField.__init__( self, tuple(valence), a.dim, g=g )
        self.a = a
        self.pairs = pairs
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        return self.a(*xs).contract( *self.pairs )

class TransposeTensorField(TensorField):
    def __init__( self, tfield, perm, g=None ):
        self.tfield = tfield
        self.perm = perm
        valence = tuple(tfield.valence[i] for i in perm)
        TensorField.__init__( self, valence, tfield.dim, g=g )
    def __call__( self, *xs ):
        t = self.tfield(*xs).transpose(*self.perm)
        assert t.valence == self.valence, (t.valence, self.valence)
        return t

class TransformTensorField(TensorField):
    def __init__( self, tfield, coord_transform, g=None ):
        TensorField.__init__( self, tfield.valence, tfield.dim, g=g )
        self.tfield = tfield
        self.coord_transform = coord_transform
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        _xs = self.coord_transform(*xs)
        tensor = self.tfield(*_xs)
        if Tensor.up in self.valence:
            inverse_partial = self.coord_transform.inverse.partial(*_xs) 
        if Tensor.dn in self.valence:
            partial = self.coord_transform.partial(*xs)
        for i, updn in enumerate(self.valence):
            if updn == Tensor.up:
#                print tensor.valence, self.coord_transform.inverse.partial(*_xs).valence, "mul", (i,1)
                tensor = tensor.mul( inverse_partial, (i,1) )
#                print tensor.valence
                idxs = range(self.rank)
#                print idxs[:i],[idxs[-1]],idxs[i:-1]
                idxs = idxs[:i]+[idxs[-1]]+idxs[i:-1]
                tensor = tensor.transpose(*idxs) 
#                print tensor.valence
                assert self.valence==tensor.valence
            elif updn == Tensor.dn:
#                assert tensor.valence[i] == Tensor.dn
                tensor = tensor.mul( partial, (i,0) )
                idxs = range(self.rank)
                idxs = idxs[:i]+[idxs[-1]]+idxs[i:-1]
                tensor = tensor.transpose(*idxs) 
                assert self.valence==tensor.valence
            else:
                assert 0, "bad valence"
        return tensor

class ViewTensorField(TensorField):
    def __init__( self, tfield, coord_transform, g=None ):
        TensorField.__init__( self, tfield.valence, tfield.dim, g=g )
        self.tfield = tfield
        self.coord_transform = coord_transform
    def __call__( self, *xs ):
        assert len(xs)==self.dim
        xs = self.coord_transform(*xs)
        tensor = self.tfield(*xs)
        return tensor

class CommaTensorField(TensorField):
    def __init__( self, tfield, g=None ):
        valence = tfield.valence + (dn,)
        TensorField.__init__( self, valence, tfield.dim, g=g )
        self.tfield = tfield
    def __call__( self, *xs ):
        dim = len(xs)
        assert dim==self.dim
        if not isinstance(xs[0],MJet):
            # must be scalar .. ie. a rank-0 Jet
            xs = [ MJet(0).promote(_x) for _x in xs ]
        # now we need dim extra ranks
        rank0 = xs[0].rank
        rank = rank0 + dim
        xx = []
        for i,_x in enumerate(xs):
            _x = MJet(rank).promote(_x)
            idxs = [0]*rank
            idxs[rank0 + i] = 1
            _x[ tuple(idxs) ] = self.scalar_one
            xx.append(_x)
    
        #############################
        tensor = self.tfield( *xx ) # 
        #############################
    
        comma = Tensor( self.valence, self.dim )

        for _idxs in tensor.genidx():
            g = tensor[_idxs]
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
            comma[_idxs] = gs
#            print "comma", idxs, gs
        return comma
            

class ScalarField(ConcreteTensorField):
    def __init__( self, dim, fn, g=None ):
        ConcreteTensorField.__init__( self, (), dim, g, {():fn} )
    
class Transform(ConcreteTensorField):
    """Coordinate Transform: a 'vector' of ScalarFields """
    def __init__( self, fns, inverse=None, g=None ):
#        self.fns = [ ScalarField.promote(fn) for fn in fns ]
        dim = len(fns)
        self.dim = dim
        self.inverse = inverse
        if inverse is not None:
            assert dim == inverse.dim
            self.inverse.inverse = self
        elems = [ ((i,),fn) for i,fn in enumerate(fns) ]
        ConcreteTensorField.__init__( self, (Tensor.up,), self.dim, g, elems )
        self.partial = self.comma()

class IdentityTransform(Transform):
    def __init__( self, dim ):
        Transform.__init__( [lambda x:x]*dim, self )

