#!/usr/bin/env python

#    tensor.py : multi-dimension tensor objects.
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
from pydx.mjet import MJet, cross

up = 1; dn = -1;

def genidx( shape ):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx( shape[1:] ):
                yield (idx,)+_idx

class Tensor(object):
    """
        NB: not necessarily with the correct transformative properties.
    """
    metric = None # class attr ???
    scalar_zero = 0.0
    scalar_one = 1.0
    scalar_promote = float
    scalar_type = float
    up = 1; dn = -1;
    def __init__( self, valence = (), dim = 4, elems=None, diag=None ):
        # for now, all Tensors are "square"
        self.valence = tuple(valence) # sequence of up/dn
        for updn in self.valence:
            assert updn in (self.up,self.dn), updn
        self.rank = len(valence)
        self.shape = (dim,)*self.rank
        self.dim = dim
        self.elems = {} # map tuple to value
        # make it dense:
        if elems is None:
            for idxs in self.genidx():
                self.elems[idxs] = self.scalar_zero
        elif type(elems)==dict:
            self.elems = dict(elems)
        else:
            for idxs in self.genidx():
                elem = elems
                for idx in idxs:
                    try:
                        elem=elem[idx]
                    except:
                        raise ValueError, "bad elems"
                self.elems[idxs] = elem
        if diag is not None:
            assert len(diag) == self.dim
            for i in range(len(diag)):
                self[ (i,)*self.rank ] = diag[i]
        self.v_cache = {}
        self.v_cache[ self.valence ] = self
    def identity( cls, valence = (), dim = 4 ):
        tensor = cls( valence, dim )
        for idxs in tensor.genidx():
            on_diag = True
            if idxs:
                i = idxs[0]
            for j in idxs:
                if j!=i:
                    on_diag = False
            if on_diag:
                tensor[idxs] = cls.scalar_one
            else:
                tensor[idxs] = cls.scalar_zero
        return tensor
#        diag = [self.scalar_one for i in range(len(valence))]
#        tensor = cls( valence, dim, diag )
#        return tensor
    identity = classmethod(identity)
    def zero( self ):
        tensor = self.__class__( self.valence, dim=self.dim )
        return tensor
    def get_valence( self, v ):
        """
            raise or lower indices (get a new valence).
        """
        for updn in v:
            assert updn in (self.up,self.dn), "bad valence: %s"%updn
        assert len(v)==self.rank
        return self.v_cache[v]
    def __getitem__( self, idxs ):
#        print self.elems, "__getitem__", idxs
        if type(idxs)==int:
            idxs = idxs,
        if type(idxs)==slice:
            start, stop, step = idxs.indices(len(self))
            idxs = range(start,stop,step)
            tensor = Tensor(self.valence, len(idxs))
#            print "start, stop, step", start, stop, step  
            for idx in idxs:
                tensor[idx-start] = self[idx]
            return tensor
        else:
            assert type(idxs)==tuple, "what's this: %s"%repr(idxs)
            assert len(idxs)<=self.rank
            assert len(idxs)==self.rank, "Not implemented"
            value = self.elems.get( idxs, None )
            if value is None:
                raise IndexError, idxs
            return value
#        return self.elems.get( idxs, self.scalar_zero )
    def __len__( self ):
        return self.dim
    def __setitem__( self, idxs, val ):
        assert val is not None
#        print "__setitem__", idxs, val
        if type(idxs)==int:
            idxs = idxs,
#        if idxs==slice(None):
#            idxs = tuple(range(len(val)))
#        if type(idxs)==slice:
#            tensor = Tensor(self.valence, self.dim)
#            start, stop, step = idxs.indices(len(self))
#            for idx in range(start,stop,step):
#                tensor[idx-start] = self[idx]
#            return tensor
#        else:
        assert type(idxs)==tuple
        free_rank = self.rank-len(idxs)
        assert free_rank >= 0, "set valence %s at %s "%( self.valence, idxs )
        if free_rank == 0:
            self.elems[idxs] = val
        elif free_rank == 1:
            assert len(val)==self.dim
            for i in range(self.dim):
                self.elems[idxs+(i,)] = val[i]
        else:
            for _idxs in cross((self.dim)*free_rank):
                self.elems[idxs+_idxs] = val[_idxs]
    def scalar( self ):
        tensor = Tensor( self.valence, self.dim )
        for idxs in self.genidx():
            r = self[idxs]
            r = MJet(0).promote(r)[()]
            tensor[idxs] = r
        return tensor
    def is_close( self, other, epsilon=1e-10 ):
        assert self.valence==other.valence
        assert self.dim==other.dim
        for idxs in self.genidx():
            r = self[idxs]-other[idxs]
            r = MJet(0).promote(r)[()]
            if abs(r)>epsilon:
                return False
        return True
    def idxstr( self, idxs, onebased=0 ):
        _v = None
        components = []
        for i, v in enumerate(self.valence):
            s = str(idxs[i]+onebased)
            if _v is not v:
                s = {up:"^",dn:"_"}[v]+s
                _v = v
            components.append( s )
        return ''.join( components ) 
    def str( self, name="" ):
        lines = []
        for idxs in self.genidx():
            lines.append( "%s%s %s" % ( name, self.idxstr(idxs), str(self[idxs])) )
        return '\n'.join(lines)
    def __str__(self):
        if self.rank==0:
            return "(%s)"%(self[()])
        ss = []
        for idxs in cross( (self.dim,)*(self.rank-1) ):
            ss.append( ' '.join([ str(self[idxs+(i,)]) for i in range(self.dim) ]) )
        return '\n'.join(ss)
    __repr__ = __str__
    def clone(self):
        " deepcopy "
        tensor = self.__class__( self.valence, self.dim )
        keys = self.keys()
#        assert len(keys) == len(self.values())
        for i, value in enumerate(self.values()):
#            print keys[i], value
            tensor[keys[i]] = value
#        assert tensor.keys() == self.keys(), tensor.keys()
        return tensor
    def __contains__( self, key ):
        return key in self.elems
    def keys( self ):
        keys = self.elems.keys()
        keys.sort()
        return keys
    def genidx( self ):
        " all possible keys "
        if len(self.shape)==0:
            yield ()
        else:
            for idx in range(self.shape[0]):
                for _idx in genidx( self.shape[1:] ):
                    yield (idx,)+_idx
    def __add__( self, other ):
        tensor = self.zero()
        for idx in self.genidx():
            tensor[idx] = self[idx]+other[idx]
        return tensor
    def __iadd__( self, other ):
        for idx in self.genidx():
            self[idx] += other[idx]
        return self
    def __sub__( self, other ):
        tensor = self.zero()
        for idx in self.genidx():
            tensor[idx] = self[idx]-other[idx]
        return tensor
    def __isub__( self, other ):
        for idx in self.genidx():
            self[idx] -= other[idx]
        return self
    def __rmul__( self, other ):
        tensor = self.zero()
        for idx in self.genidx():
            tensor[idx] = other * self[idx]
        return tensor
    def __radd__( self, other ):
        tensor = self.zero()
        for idx in self.genidx():
            tensor[idx] = other + self[idx]
        return tensor

    def contract( self, *pairs ):
        " pairs: sequence of self index, self index "
        for i, j in pairs:
            assert self.valence[i] != self.valence[j] # pair up with dn, and vice-versa
            assert i not in [j for i,j in pairs]
            assert j not in [i for i,j in pairs]
#        src_idx = {}; tgt_idx = {}; # map indexes
        valence = []
        left = [i for i,j in pairs]
        right = [j for i,j in pairs]
        contracted = left+right # the contracted indexs
        for idx, v in enumerate(self.valence):
            if idx not in contracted:
                # we are not contracting this idx
                valence.append( v )
#                src_idx[idx] = len(valence)-1 # this index in self maps to this index in tensor
        tensor = self.__class__( valence, self.dim )
        assert ( self.rank - tensor.rank ) % 2 == 0, "internal error"
        kernel = (self.dim,) * (( self.rank - tensor.rank )/2)
        assert len(kernel) == len(pairs)
#        print "contracted:", contracted
        for tgt_idx in tensor.genidx():
            res = self.scalar_zero
#            print "tgt_idx:", tgt_idx
            for kernel_idx in genidx( kernel ):
#                print "  kernel_idx:", kernel_idx
                # weave _idx into tgt_idx to get src_idx
                src_idx = []
                i = 0 # index into tgt_idx
                for srci in range(self.rank):
#                    print "    srci:", srci
                    if srci in left:
                        kernel_i = left.index(srci)
                        src_idx.append( kernel_idx[kernel_i] )
                    elif srci in right:
                        kernel_i = right.index(srci)
                        src_idx.append( kernel_idx[kernel_i] )
                    else:
                        src_idx.append( tgt_idx[i] )
                        i += 1
                assert i == len(tgt_idx)
                assert len(src_idx) == self.rank
                src_idx = tuple(src_idx)
                if src_idx in self:
                    res = res + self[src_idx]
#            if res != self.scalar_zero:
            tensor[ tgt_idx ] = res
        return tensor

    def outer( self, other ):
        assert self.dim == other.dim
        valence = self.valence + other.valence
        tensor = self.__class__( valence, self.dim )
        for idx in self.genidx():
            for _idx in other.genidx():
#                if self[idx] != self.scalar_zero and other[_idx] != self.scalar_zero:
                tensor[ idx + _idx ] = self[idx]*other[_idx]
        return tensor
    __mul__ = outer

    def mul( self, other, *pairs ):
        " pairs: sequence of self index, other index "
        for i, j in pairs:
            # pair up with dn, and vice-versa
            assert self.valence[i] != other.valence[j], (self.valence, i, other.valence, j)
        assert self.dim == other.dim
        tensor = self.outer(other) # big
        pairs = [ (i, j+self.rank) for i,j in pairs ]
        tensor = tensor.contract( *pairs )
        return tensor

    def transpose( self, *axes ):
        assert len(axes)==self.rank
        valence = [self.valence[i] for i in axes]
        tensor = Tensor(valence, self.dim)
        for idx in self.genidx():
            _idx = tuple([idx[axes[i]] for i in range(self.rank)])
            tensor[idx] = self[_idx] #.clone()
        return tensor

    def transform( self, transform=None, inverse=None ):
        tensor = self
        for i, updn in enumerate(self.valence):
            if updn == Tensor.up:
                tensor = tensor.mul( inverse, (i,1) )
            elif updn == Tensor.dn:
                tensor = tensor.mul( transform, (i,0) )
            else:
                assert 0, "bad valence"
        return tensor

    def apply( self, func ):
        for idx in self.genidx():
            self[idx] = func( self[idx] )
        return self
scalar.Scalar.clients.append(Tensor)

def Scalar( value ):
    tensor = Tensor()
    tensor[()] = value
    return tensor

#
###################################
#
