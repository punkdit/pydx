#!/usr/bin/env python

#    mpfi.pyx: wrappers for mpfi lib
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



cdef extern from *:
    ctypedef unsigned int __c_size_t "size_t"

cdef extern from "errno.h":
    int errno
    void perror(char *s)
    char *strerror(int errnum)

cdef extern from "stdio.h":
    ctypedef void FILE
    FILE *stdout
    FILE *stderr

cdef extern from "string.h":
    void *__c_memcpy "memcpy" ( void*, void*, __c_size_t )
    int __c_memcmp "memcmp" ( void*, void*, __c_size_t )

cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void free(void*)


cdef extern from "Python.h":
    void*  PyMem_Malloc(__c_size_t)
    void*  PyMem_Realloc(void *p, __c_size_t n)
    void   PyMem_Free(void *p)

cdef extern from "gmp.h":
    ctypedef int mp_prec_t
    ctypedef int mp_exp_t
    ctypedef unsigned int mp_limb_t
    struct __mpf_struct:
        int _mp_prec  #/* Max precision, in number of `mp_limb_t's.  Set by mpf_init and modified by mpf_set_prec.  The area pointed to by the _mp_d field contains `prec' + 1 limbs.  */
        int _mp_size  # /* abs(_mp_size) is the number of limbs the last field points to.  If _mp_size is negative this is a negative number.  */
        mp_exp_t _mp_exp  #/* Exponent, in the base of `mp_limb_t'.  */
        mp_limb_t *_mp_d  #/* Pointer to the limbs.  */

    ctypedef __mpf_struct *mpf_t
    ctypedef __mpf_struct *mpf_ptr
    ctypedef __mpf_struct *mpf_srcptr
    void mpf_init( mpf_ptr )
    void mpf_set( mpf_t, mpf_t )
    void mpf_set_d( mpf_t, double )
    int mpf_get_prec( mpf_t )
    void __c_mpf_set_default_prec "mpf_set_default_prec"( int )
    double mpf_get_d( mpf_srcptr )
    __c_size_t mpf_out_str (FILE *, int, __c_size_t, mpf_srcptr)
    char *mpf_get_str  (char *, mp_exp_t *, int, size_t, mpf_srcptr)


cdef extern from "gmpy.h":
    void import_gmpy()

    ctypedef struct PympfObject:
        mpf_t f
        int rebits
    object Pympf_new(int bits)
    int Pympf_Check(o)


cdef extern from "mpfr.h":
    ctypedef enum mpfr_rnd_t:
        GMP_RNDN # nearest
        GMP_RNDZ # zero
        GMP_RNDU # up
        GMP_RNDD # down
        GMP_RND_MAX
        GMP_RNDNA

    ctypedef int mpfr_prec_t
    ctypedef int mpfr_sign_t

    ctypedef struct __mpfr_struct:
        mpfr_prec_t  _mpfr_prec
        mpfr_sign_t  _mpfr_sign
        mp_exp_t     _mpfr_exp
        mp_limb_t   *_mpfr_d
    
    ctypedef __mpfr_struct mpfr_t[1]
#    ctypedef __mpfr_struct *mpfr_t
    ctypedef __mpfr_struct *mpfr_ptr
    ctypedef __mpfr_struct *mpfr_srcptr


    void mpfr_set_default_prec(mp_prec_t)
    mp_prec_t mpfr_get_default_prec()

    void mpfr_init( mpfr_ptr )
    void mpfr_init2( mpfr_ptr, mpfr_prec_t )
    void mpfr_clear( mpfr_ptr )

    mpfr_prec_t mpfr_get_prec( mpfr_srcptr )

    int mpfr_set_f (mpfr_ptr, mpf_srcptr, mpfr_rnd_t)
    int mpfr_get_f(mpf_ptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_set_d (mpfr_ptr, double, mpfr_rnd_t)
    double mpfr_get_d (mpfr_srcptr, mpfr_rnd_t)

    int mpfr_pow (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_pow_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
    int mpfr_add (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_sub (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_mul (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_div (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)

    char*mpfr_get_str (char*, mp_exp_t*, int, __c_size_t, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_strtofr (mpfr_t rop, char *nptr, char **endptr, int base, mpfr_rnd_t rnd)
    void mpfr_free_str (char *str)

    int mpfr_check_range (mpfr_t X, int T, mpfr_rnd_t RND)

    int mpfr_underflow_p ()
    int mpfr_overflow_p ()
    int mpfr_nanflag_p ()
    int mpfr_inexflag_p ()
    int mpfr_erangeflag_p ()

    void mpfr_clear_flags ()



cdef extern from "mpfi.h":
    
    ctypedef struct __mpfi_struct:
        __mpfr_struct left
        __mpfr_struct right
    
    ctypedef __mpfi_struct mpfi_t[1]
#    ctypedef __mpfi_struct *mpfi_t
    ctypedef __mpfi_struct *mpfi_ptr
    ctypedef __mpfi_struct *mpfi_srcptr
    
    
    # Rounding                                     
    int mpfi_round_prec(mpfi_ptr,mp_prec_t prec)
    
    
    # Initialization, destruction and assignment   
    
    # initializations 
    void   mpfi_init       (mpfi_ptr)
    void   mpfi_init2      (mpfi_ptr, mp_prec_t)
    
    void   mpfi_clear      (mpfi_ptr)
    
    # mpfi bounds have the same precision 
    mp_prec_t mpfi_get_prec(mpfi_srcptr)
    void mpfi_set_prec(mpfi_ptr,mp_prec_t)
    
    
    # assignment functions                         
    int   mpfi_set        (mpfi_ptr, mpfi_srcptr)
    int   mpfi_set_si     (mpfi_ptr, long)
    int   mpfi_set_ui     (mpfi_ptr, unsigned long)
    int   mpfi_set_d      (mpfi_ptr, double)
    int   mpfi_set_z      (mpfi_ptr, mpz_srcptr)
    int   mpfi_set_q      (mpfi_ptr, mpq_srcptr)
    int   mpfi_set_fr   (mpfi_ptr, mpfr_srcptr)
    int   mpfi_set_str    (mpfi_ptr, char *, int)
    
    # combined initialization and assignment functions 
    int   mpfi_init_set        (mpfi_ptr, mpfi_srcptr)
    int   mpfi_init_set_si     (mpfi_ptr, long)
    int   mpfi_init_set_ui     (mpfi_ptr, unsigned long)
    int   mpfi_init_set_d      (mpfi_ptr, double)
    int   mpfi_init_set_z      (mpfi_ptr, mpz_srcptr)
    int   mpfi_init_set_q      (mpfi_ptr, mpq_srcptr)
    int   mpfi_init_set_fr   (mpfi_ptr, mpfr_srcptr)
    int   mpfi_init_set_str    (mpfi_ptr, char *, int)
    
    # swapping two intervals 
    void mpfi_swap (mpfi_ptr, mpfi_ptr)
    
    
    # Various useful interval functions            
    # with scalar or interval results              
    
    # absolute diameter                            
    int mpfi_diam_abs(mpfr_ptr, mpfi_srcptr)
    # relative diameter                            
    int mpfi_diam_rel(mpfr_ptr, mpfi_srcptr)
    # diameter: relative if the interval does not contain 0 
    # absolute otherwise                                    
    int mpfi_diam(mpfr_ptr, mpfi_srcptr)
    # magnitude: the largest absolute value of any element 
    int mpfi_mag(mpfr_ptr, mpfi_srcptr)
    # mignitude: the smallest absolute value of any element 
    int mpfi_mig(mpfr_ptr, mpfi_srcptr)
    # middle of y                                           
    int mpfi_mid (mpfr_ptr, mpfi_srcptr)
    # picks randomly a point m in y 
    void mpfi_alea (mpfr_ptr, mpfi_srcptr)
    
    # Conversions                                  
    double mpfi_get_d (mpfi_srcptr)
    void mpfi_get_fr (mpfr_ptr,mpfi_srcptr)
    
    # Basic arithmetic operations                  
    
    # arithmetic operations between two interval operands 
    int   mpfi_add        (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int   mpfi_sub        (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int   mpfi_mul        (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int   mpfi_div        (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and a double prec. floating-point 
    int   mpfi_add_d      (mpfi_ptr, mpfi_srcptr, double)
    int   mpfi_sub_d      (mpfi_ptr, mpfi_srcptr, double)
    int   mpfi_d_sub      (mpfi_ptr, double, mpfi_srcptr)
    int   mpfi_mul_d      (mpfi_ptr, mpfi_srcptr, double)
    int   mpfi_div_d      (mpfi_ptr, mpfi_srcptr, double)
    int   mpfi_d_div      (mpfi_ptr, double, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and an unsigned long integer 
    int   mpfi_add_ui     (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_sub_ui     (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_ui_sub     (mpfi_ptr, unsigned long, mpfi_srcptr)
    int   mpfi_mul_ui     (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_div_ui     (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_ui_div     (mpfi_ptr, unsigned long, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and a long integer 
    int   mpfi_add_si     (mpfi_ptr, mpfi_srcptr, long)
    int   mpfi_sub_si     (mpfi_ptr, mpfi_srcptr, long)
    int   mpfi_si_sub     (mpfi_ptr, long, mpfi_srcptr)
    int   mpfi_mul_si     (mpfi_ptr, mpfi_srcptr, long)
    int   mpfi_div_si     (mpfi_ptr, mpfi_srcptr, long)
    int   mpfi_si_div     (mpfi_ptr, long, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and a multiple prec. integer 
    int   mpfi_add_z     (mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int   mpfi_sub_z     (mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int   mpfi_z_sub     (mpfi_ptr, mpz_srcptr, mpfi_srcptr)
    int   mpfi_mul_z     (mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int   mpfi_div_z     (mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int   mpfi_z_div     (mpfi_ptr, mpz_srcptr, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and a multiple prec. rational 
    int   mpfi_add_q     (mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int   mpfi_sub_q     (mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int   mpfi_q_sub     (mpfi_ptr, mpq_srcptr, mpfi_srcptr)
    int   mpfi_mul_q     (mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int   mpfi_div_q     (mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int   mpfi_q_div     (mpfi_ptr, mpq_srcptr, mpfi_srcptr)
    
    # arithmetic operations between an interval operand and a mult. prec. floating-pt nb 
    int   mpfi_add_fr  (mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int   mpfi_sub_fr  (mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int   mpfi_fr_sub  (mpfi_ptr, mpfr_srcptr, mpfi_srcptr)
    int   mpfi_mul_fr  (mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int   mpfi_div_fr  (mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int   mpfi_fr_div  (mpfi_ptr, mpfr_srcptr, mpfi_srcptr)
    
    # arithmetic operations taking a single interval operand 
    int   mpfi_neg        (mpfi_ptr, mpfi_srcptr)
    int   mpfi_sqr        (mpfi_ptr ,mpfi_srcptr)
    # the inv function generates the whole real interval if 0 is in the interval defining the divisor 
    int   mpfi_inv        (mpfi_ptr, mpfi_srcptr)
    # the sqrt of a (partially) negative interval is a NaN 
    int   mpfi_sqrt       (mpfi_ptr, mpfi_srcptr)
    # the first interval contains the absolute values of 
    # every element of the second interval 
    int   mpfi_abs        (mpfi_ptr, mpfi_srcptr)
    
    # various operations 
    int   mpfi_mul_2exp   (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_mul_2ui    (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_mul_2si    (mpfi_ptr, mpfi_srcptr, long)
    int   mpfi_div_2exp   (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_div_2ui    (mpfi_ptr, mpfi_srcptr, unsigned long)
    int   mpfi_div_2si    (mpfi_ptr, mpfi_srcptr, long)
    
    # Special functions                                        
    int mpfi_log  (mpfi_ptr, mpfi_srcptr)
    int mpfi_exp  (mpfi_ptr, mpfi_srcptr)
    int mpfi_exp2 (mpfi_ptr, mpfi_srcptr)
    
    int mpfi_cos  (mpfi_ptr, mpfi_srcptr)
    int mpfi_sin  (mpfi_ptr, mpfi_srcptr)
    int mpfi_tan  (mpfi_ptr, mpfi_srcptr)
    int mpfi_acos (mpfi_ptr, mpfi_srcptr)
    int mpfi_asin (mpfi_ptr, mpfi_srcptr)
    int mpfi_atan (mpfi_ptr, mpfi_srcptr)
    
    int mpfi_cosh (mpfi_ptr, mpfi_srcptr)
    int mpfi_sinh (mpfi_ptr, mpfi_srcptr)
    int mpfi_tanh (mpfi_ptr, mpfi_srcptr)
    int mpfi_acosh (mpfi_ptr, mpfi_srcptr)
    int mpfi_asinh (mpfi_ptr, mpfi_srcptr)
    int mpfi_atanh (mpfi_ptr, mpfi_srcptr)
    
    int mpfi_log1p (mpfi_ptr, mpfi_srcptr)
    int mpfi_expm1 (mpfi_ptr, mpfi_srcptr)
    
    int mpfi_log2 (mpfi_ptr, mpfi_srcptr)
    int mpfi_log10 (mpfi_ptr, mpfi_srcptr)
    
    int mpfi_const_log2(mpfi_ptr)
    int mpfi_const_pi(mpfi_ptr)
    int mpfi_const_euler(mpfi_ptr)
    
    # Comparison functions                                     
    # Warning: the meaning of interval comparison is not clearly defined 
    # customizable comparison functions 
    
    int    (*mpfi_cmp)     (mpfi_srcptr,mpfi_srcptr)
    int    (*mpfi_cmp_d)   (mpfi_srcptr,double)
    int    (*mpfi_cmp_ui)  (mpfi_srcptr,unsigned long)
    int    (*mpfi_cmp_si)  (mpfi_srcptr,long)
    int    (*mpfi_cmp_z)   (mpfi_srcptr,mpz_srcptr)
    int    (*mpfi_cmp_q)   (mpfi_srcptr,mpq_srcptr)
    int    (*mpfi_cmp_fr)(mpfi_srcptr,mpfr_srcptr)
    
    int    (*mpfi_is_pos)     (mpfi_srcptr)
    int    (*mpfi_is_nonneg)  (mpfi_srcptr)
    int    (*mpfi_is_neg)     (mpfi_srcptr)
    int    (*mpfi_is_nonpos)  (mpfi_srcptr)
    int    (*mpfi_is_zero)    (mpfi_srcptr)
    int    (*mpfi_is_strictly_pos) (mpfi_srcptr)
    int    (*mpfi_is_strictly_neg) (mpfi_srcptr)
    
    # default comparison functions 
    int    mpfi_is_pos_default          (mpfi_srcptr)
    int    mpfi_is_nonneg_default       (mpfi_srcptr)
    int    mpfi_is_neg_default          (mpfi_srcptr)
    int    mpfi_is_nonpos_default       (mpfi_srcptr)
    int    mpfi_is_zero_default         (mpfi_srcptr)
    int    mpfi_is_strictly_neg_default (mpfi_srcptr a)
    int    mpfi_is_strictly_pos_default (mpfi_srcptr a)
    
    int    mpfi_cmp_default      (mpfi_srcptr,mpfi_srcptr)
    int    mpfi_cmp_d_default    (mpfi_srcptr,double)
    int    mpfi_cmp_ui_default   (mpfi_srcptr,unsigned long)
    int    mpfi_cmp_si_default   (mpfi_srcptr,long)
    int    mpfi_cmp_z_default    (mpfi_srcptr,mpz_srcptr)
    int    mpfi_cmp_q_default    (mpfi_srcptr,mpq_srcptr)
    int    mpfi_cmp_fr_default (mpfi_srcptr,mpfr_srcptr)
    
    
    int mpfi_has_zero (mpfi_srcptr)
    
    int mpfi_nan_p (mpfi_srcptr)
    int mpfi_inf_p (mpfi_srcptr)
    int mpfi_bounded_p (mpfi_srcptr)
    
    # Interval manipulation 
    
    # operations related to the internal representation by endpoints 
    
    # get left or right bound of the interval defined by the second argument and put the result in the first one 
    int   mpfi_get_left   (mpfr_ptr,mpfi_srcptr)
    int   mpfi_get_right  (mpfr_ptr,mpfi_srcptr)
    
    int   mpfi_revert_if_needed  (mpfi_ptr)
    
    # Set operations on intervals 
    # "Convex hulls" 
    # extends the interval defined by the first argument so that it contains the second one 
    
    int   mpfi_put        (mpfi_ptr,mpfi_srcptr)
    int   mpfi_put_d      (mpfi_ptr,double)
    int   mpfi_put_si     (mpfi_ptr,long)
    int   mpfi_put_ui     (mpfi_ptr,unsigned long)
    int   mpfi_put_z      (mpfi_ptr,mpz_srcptr)
    int   mpfi_put_q      (mpfi_ptr,mpq_srcptr)
    int   mpfi_put_fr     (mpfi_ptr,mpfr_srcptr)
    
    # builds an interval whose left bound is the lower (round -infty) than the second argument and the right bound is greater (round +infty) than the third one 
    
    int   mpfi_interv_d   (mpfi_ptr,double,double)
    int   mpfi_interv_si  (mpfi_ptr,long,long)
    int   mpfi_interv_ui  (mpfi_ptr,unsigned long,unsigned long)
    int   mpfi_interv_z   (mpfi_ptr,mpz_srcptr,mpz_srcptr)
    int   mpfi_interv_q   (mpfi_ptr,mpq_srcptr,mpq_srcptr)
    int   mpfi_interv_fr  (mpfi_ptr,mpfr_srcptr,mpfr_srcptr)
    
    # Inclusion tests 
    # tests if the first argument is inside the interval defined by the second one 
    int   mpfi_is_strictly_inside   (mpfi_srcptr,mpfi_srcptr)
    int   mpfi_is_inside        	(mpfi_srcptr, mpfi_srcptr)
    int   mpfi_is_inside_d      	(double, mpfi_srcptr)
    int   mpfi_is_inside_ui     	(unsigned long, mpfi_srcptr)
    int   mpfi_is_inside_si     	(long, mpfi_srcptr)
    int   mpfi_is_inside_z      	(mpz_srcptr,mpfi_srcptr)
    int   mpfi_is_inside_q      	(mpq_srcptr,mpfi_srcptr)
    int   mpfi_is_inside_fr   	(mpfr_srcptr,mpfi_srcptr)
    
    # set operations 
    int   mpfi_is_empty         (mpfi_srcptr)
    int   mpfi_intersect        (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int   mpfi_union            (mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    
    # complement... : to do later 
    
    
    # Miscellaneous 
    
    # adds the second argument to the right bound of the first one and subtracts the second argument to the left bound of the first one 
    int  mpfi_increase         (mpfi_ptr,mpfr_srcptr)
    # keeps the same center and multiply the radius by 2*(1+fact) 
    int mpfi_blow(mpfi_ptr, mpfi_srcptr, double)
    # splits the interval into 2 halves 
    int mpfi_bisect(mpfi_ptr, mpfi_ptr, mpfi_srcptr)
    
    char * mpfi_get_version()
    
    # Error handling 
    
    void   mpfi_reset_error()
    void   mpfi_set_error(int)
    int    mpfi_is_error()
    
    #define MPFI_FLAGS_BOTH_ENDPOINTS_EXACT       0
    #define MPFI_FLAGS_LEFT_ENDPOINT_INEXACT      1
    #define MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT     2
    #define MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT     3
    
    #define MPFI_BOTH_ARE_EXACT(x) ( (int)(x) == 0 )
    #define MPFI_LEFT_IS_INEXACT(x) ( (int)(x)%2 )
    #define MPFI_RIGHT_IS_INEXACT(x) ( (int)(x)/2 )
    #define MPFI_BOTH_ARE_INEXACT(x) ( (int)(x) == 3 )
    
    #define MPFI_REVERT_INEXACT_FLAGS(x) ( ((x)==1) ? 2 : ((x)==2) ? 1 : x )
    
    #define MPFI_NAN_P(a) ( MPFR_IS_NAN(&(a->left)) || MPFR_IS_NAN (&(a->right)) )
    #define MPFI_INF_P(a) ( MPFR_IS_INF(&(a->left)) || MPFR_IS_INF (&(a->right)) )
    #define MPFI_IS_ZERO(a)  (MPFI_NAN_P(a) ? 0 : ((MPFR_SIGN(&(a->right))==0) && (MPFR_SIGN(&(a->left))==0)))
    
    #define MPFI_CLEAR(m) {mpfr_clear(&(m->right)); mpfr_clear(&(m->left));}

import_gmpy()

sizeof_mpfi = sizeof(mpfi_t)

def mpf_set_default_prec(prec):
    __c_mpf_set_default_prec(prec)

cdef set_mpf( mpf_t _x, x ):
    " set _x from x, steals a ref "
    assert Pympf_Check(x)
    mpf_set( _x, (<PympfObject*>x)[0].f )

cdef object mpf2py( mpf_t x ):
    " new ref "
    cdef object _x
    cdef int bits
    bits = mpf_get_prec( x )
    _x = Pympf_new( bits )
    mpf_set( (<PympfObject*>_x)[0].f, x )
    return <object>_x

cdef object mpfr2py( mpfr_t x, mpfr_rnd_t rnd ):
    " new ref "
    cdef object _x
    cdef int bits
    bits = mpfr_get_prec( x )
    _x = Pympf_new( bits )
    mpfr_get_f( (<PympfObject*>_x)[0].f, x, rnd )
    return <object>_x

import gmpy
from gmpy import mpf

cdef object mpf_type
mpf_type = type(mpf(1.0))

cdef class Interval

cdef object c_promote(x):
    if isinstance(x,Interval):
        pass
    elif hasattr(x,"__float__") or Pympf_Check(x):
        x = Interval(x)
    else:
        x = None
    return x

def promote(x):
    if isinstance(x,Interval):
        return x
    if hasattr(x,"__float__") or Pympf_Check(x):
        return Interval(x)
    return None

def set_default_prec(prec):
    mpfr_set_default_prec(prec)
def get_default_prec():
    return mpfr_get_default_prec()

class DisorderException(Exception):
    pass

cdef class Interval:
    cdef mpfi_t x
    def __init__( self, lower=None, upper=None ):
        cdef Interval _lower
        cdef mpfr_t fr_lower
        cdef mpfr_t fr_upper
        cdef PympfObject *pympf_lower
        cdef PympfObject *pympf_upper

        if isinstance(lower, Interval):
            _lower = lower
            mpfi_init_set( self.x, _lower.x )
            assert upper is None, "but first arg is already an Interval ?!?!"
        else:
            if hasattr(lower,"__float__"):
                mpfi_init_set_d( self.x, lower )
            elif Pympf_Check(lower):
                pympf_lower = <PympfObject*>lower
                mpfr_init(fr_lower)
                mpfr_set_f(fr_lower,pympf_lower[0].f,GMP_RNDD)
                mpfi_init_set_fr( self.x, fr_lower )
                mpfr_clear(fr_lower)
            elif lower is None:
                mpfi_init( self.x )
                assert upper is None
            else:
                mpfi_init( self.x )
                raise TypeError, "what's this: %s"%repr(lower)
            if hasattr(upper,"__float__"):
                mpfi_put_d( self.x, upper )
            elif Pympf_Check(upper):
                pympf_upper = <PympfObject*>upper
                mpfr_init(fr_upper)
                mpfr_set_f(fr_upper,pympf_upper[0].f,GMP_RNDU)
                mpfi_put_fr( self.x, fr_upper )
                mpfr_clear(fr_upper)
            elif upper is not None:
                mpfi_init( self.x )
                raise TypeError, "what's this: %s"%upper
        assert (lower is None and upper is None) or not mpfi_is_empty( self.x ),\
            (self.lower.digits(),self.upper.digits(),lower,upper)
    def __dealloc__( self ):
        mpfi_clear( self.x )
    def get_addr( self ):
        cdef long addr
        addr = <long>self.x
        return addr

    property lower:
        def __get__( self ):
            cdef mpfr_t fr_lower
            mpfr_init( fr_lower )
            mpfi_get_left( fr_lower, self.x )
            lower = mpfr2py( fr_lower, GMP_RNDD )
            mpfr_clear( fr_lower )
            return lower
    property upper:
        def __get__( self ):
            cdef mpfr_t fr_upper
            mpfr_init( fr_upper )
            mpfi_get_right( fr_upper, self.x )
            upper = mpfr2py( fr_upper, GMP_RNDU )
            mpfr_clear( fr_upper )
            return upper
    property prec:
        def __get__( self ):
            return mpfi_get_prec( self.x )
        def __set__( self, prec ):
            mpfi_round_prec( self.x, prec )

    def __repr__( self ):
        return self.to_string()
#        return "Interval(%s,%s)" % (self.lower, self.upper)
    def __str__( self ):
#        return "Interval(%.6f,%.6f)" % (self.lower, self.upper)
        return "[%.9f,%.9f]" % (self.lower, self.upper)
    def to_string( self ):
        # char * mpfr_get_str (char *str, mp_exp_t *expptr, int b, size_t n, mpfr_t op, mp_rnd_t rnd)
        cdef char *_s_lower, *_s_upper
        cdef mp_exp_t exp
        cdef mpfr_t fr_lower, fr_upper

        mpfr_init( fr_lower )
        mpfi_get_left( fr_lower, self.x )
        mpfr_init( fr_upper )
        mpfi_get_right( fr_upper, self.x )

        _s_lower = mpfr_get_str(NULL, &exp, 10, 0, fr_lower, GMP_RNDD )
        s_lower = _s_lower
        if s_lower is None:
            raise Exception
        if s_lower[0]=='-':
            sl = s_lower[:2]+'.'+s_lower[2:]+'e'+str(exp-1)
        else:
            sl = s_lower[:1]+'.'+s_lower[1:]+'e'+str(exp-1)

        _s_upper = mpfr_get_str(NULL, &exp, 10, 0, fr_upper, GMP_RNDD )
        s_upper = _s_upper
        if s_upper is None:
            raise Exception
        if s_upper[0]=='-':
            su = s_upper[:2]+'.'+s_upper[2:]+'e'+str(exp-1)
        else:
            su = s_upper[:1]+'.'+s_upper[1:]+'e'+str(exp-1)

        mpfr_free_str(_s_lower)
        mpfr_free_str(_s_upper)
        mpfr_clear( fr_lower )
        mpfr_clear( fr_upper )
        return "Interval(mpf('%s'), mpf('%s'))" % (sl, su)
    def __hash__( self ):
        return hash(self.to_string())

    def __richcmp__( _x, _y, op ):
        cdef Interval x, y, z
        x = c_promote(_x)
        if x is None:
            return NotImplemented
        y = c_promote(_y)
        if y is None:
            return NotImplemented
        if op == 0: # <
            if x.upper < y.lower:
                return True
        elif op == 1: # <=
            if x.upper <= y.lower:
                return True
        elif op == 2: # ==
            if x.lower == y.lower and x.upper == y.upper:
                return True
        elif op == 3: # !=
            if x.lower != y.lower or x.upper != y.upper:
                return True
        elif op == 4: # >
            if x.lower > y.upper:
                return True
        elif op == 5: # >=
            if x.lower >= y.upper:
                return True
        return False

    def freeze(self):
        return self

    def __neg__( x ):
        cdef Interval y
        y = Interval()
        mpfi_neg( y.x, x.x )
        if mpfi_is_empty( y.x ):
            raise DisorderException
        return y
    def __pos__( x ):
        return x

    def __abs__( x ):
        cdef Interval y
        y = Interval()
        mpfi_abs( y.x, x.x )
        if mpfi_is_empty( y.x ):
            raise DisorderException
        return y

    def __add__( _x, _y ):
        cdef Interval x, y, z
        x = c_promote(_x)
        if x is None:
            return NotImplemented
        y = c_promote(_y)
        if y is None:
            return NotImplemented
        z = Interval()
        mpfi_add( z.x, x.x, y.x )
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z

    def __sub__( _x, _y ):
        cdef Interval x, y, z
        x = c_promote(_x)
        if x is None:
            return NotImplemented
        y = c_promote(_y)
        if y is None:
            return NotImplemented
        z = Interval()
        mpfi_sub( z.x, x.x, y.x )
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z
    def __mul__( _x, _y ):
        cdef Interval x, y, z
        x = c_promote(_x)
        if x is None:
            return NotImplemented
        y = c_promote(_y)
        if y is None:
            return NotImplemented
        z = Interval()
        mpfi_mul( z.x, x.x, y.x )
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z
    def __div__( _x, _y ):
        cdef Interval x, y, z
        x = c_promote(_x)
        if x is None:
            return NotImplemented
        y = c_promote(_y)
        if y is None:
            return NotImplemented
        z = Interval()
        mpfi_div( z.x, x.x, y.x )
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z
    def __pow__( _x, y, w ):
        cdef Interval x, z
        cdef mpfr_t x_lower, x_upper
        cdef mpfr_t z_lower, z_upper

        x = c_promote(_x)
        if x is None:
            return NotImplemented
        if type(y)==Interval and y.width()==0.0 and int(y.lower)==y.lower:
            doweusethisatall
            y=int(y.lower)
        if type(y) in (int,float):
            if int(y)==y:
                y = int(y)
                if y==0:
                    return Interval(1.0)
                if y==1:
                    return x
                if y==2:
                    return x.sqr()
                mpfr_init( x_lower )
                mpfr_init( x_upper )
                mpfr_init( z_lower )
                mpfr_init( z_upper )
                mpfi_get_left( x_lower, x.x )
                mpfi_get_right( x_upper, x.x )
                z = Interval()
                mpfi_get_left( z_lower, z.x )
                mpfi_get_right( z_upper, z.x )
                mpfr_pow_si( z_lower, x_lower, y, GMP_RNDD )
                mpfr_pow_si( z_upper, x_upper, y, GMP_RNDU )
                mpfi_set_fr( z.x, z_lower )
                mpfi_put_fr( z.x, z_upper )
                mpfr_clear( x_lower )
                mpfr_clear( x_upper )
                mpfr_clear( z_lower )
                mpfr_clear( z_upper )
                if mpfi_is_empty( z.x ):
                    raise DisorderException((x, y, z))
                return z
            elif y==0.5:
                return _x.sqrt() 
            z = (y*x.log()).exp()
            if mpfi_is_empty( z.x ):
                raise DisorderException((x, y, z))
            return z
        y = c_promote(y)
        if y is None:
            return NotImplemented
        assert 0, "implement me! "
        return NotImplemented
        # return Interval(x.lower**y,x.upper**y) # BROKEN
        # int mpfr_pow (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
        # this sucks too much:
        #  z = Interval(1.0)
        #  while y:
        #    z = z*x
        #    y = y - 1
        #  return z
        
    def width(x):
        cdef mpfr_t fr_y
        cdef mpf_t f_y
        cdef int bits
        cdef object py_y
        mpfr_init(fr_y)
        mpfi_diam_abs(fr_y, x.x)
        mpf_init( f_y )
        mpfr_get_f( f_y, fr_y, GMP_RNDU )
        bits = mpf_get_prec( f_y )
        py_y = Pympf_new( bits )
        mpf_set( (<PympfObject*>py_y)[0].f, f_y )
        mpfr_clear(fr_y)
        return py_y
#    def width(x):
#        return x.upper-x.lower
#    def width(x):
#        return mpfr_get_d( &x.x.right, GMP_RNDU ) - mpfr_get_d( &x.x.left, GMP_RNDD )
    def centre(x):
        cdef mpfr_t y
        mpfr_init(y)
        mpfi_mid(y, x.x)
        _y = mpfr2py(y, GMP_RNDN)
        mpfr_clear(y)
        return _y
    def overlapping(x, Interval y):
        return x.upper>=y.lower and y.upper>=x.lower
    def intersect(x, Interval y):
        cdef Interval z
        z = Interval()
        mpfi_intersect(z.x, x.x, y.x)
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z
    def hull(x, Interval y):
        cdef Interval z
        z = Interval()
        mpfi_union(z.x, x.x, y.x)
        if mpfi_is_empty( z.x ):
            raise DisorderException((x, y, z))
        return z
    def contains( x, _y ):
        cdef Interval y
        if not isinstance(_y,Interval):
            y = Interval(_y)
        else:
            y = _y
        return mpfi_is_inside( y.x, x.x )

    def is_close(self, other, EPSILON=1e-13):
        return abs(self.lower-other.lower)<EPSILON and abs(self.upper-other.upper)<EPSILON

    def bisect(x):
        cdef Interval y
        cdef Interval z
        y = Interval()
        z = Interval()
        mpfi_bisect( y.x, z.x, x.x )
        # return x,x if x width is zero ?
        if mpfi_is_empty( y.x ):
            raise DisorderException()
        if mpfi_is_empty( z.x ):
            raise DisorderException()
        return y, z

    def split(x, n):
        # subdivide x into at least n Intervals
        ivs = [x]
        while len(ivs)<n:
            _ivs = []
            for _x in ivs:
                l, r =_x.bisect()
                _ivs.append(l)
                _ivs.append(r)
            ivs = _ivs
        return ivs

#    def ceil(x):
#        return Interval(ceil(x._lower), ceil(x._upper))
#    def floor(x):
#        return Interval(floor(x._lower), floor(x._upper))
    
#    def max(x, y):
#        return Interval(max(x._lower, y._lower), max(x._upper, y._upper))
#    def min(x, y):
#        return Interval(min(x._lower, y._lower), min(x._upper, y._upper))
    
    def sqr(x):
        cdef Interval y
        y = Interval()
        mpfi_sqr( y.x, x.x )
        if mpfi_is_empty( y.x ):
            raise DisorderException()
        return y

    def sqrt( x ):
        cdef Interval y
        y = Interval()
        mpfi_sqrt( y.x, x.x )
        if mpfi_is_empty( y.x ):
            raise DisorderException()
        return y

    def hemimetric(x, y):
        ms = [
            [ abs(x.lower-y.lower), abs(x.lower-y.upper) ],
            [ abs(x.upper-y.lower), abs(x.upper-y.upper) ],
        ]
        return max(min(ms[0][0],ms[0][1]),min(ms[1][0],ms[1][1]))

    def metric( x, y ):
        " hausdorff metric "
        return max( x.hemimetric(y), y.hemimetric(x) )

    def get_pi( x ):
        cdef Interval y
        y = Interval()
        mpfi_const_pi(y.x)
        return y
    
    def exp( x ):
        cdef Interval y
        y = Interval()
        mpfi_exp( y.x, x.x )
        if mpfi_is_empty( y.x ):
            raise DisorderException()
        return y
    def log( x ):
        cdef Interval y
        y = Interval()
        mpfi_log( y.x, x.x )
        return y
    def sin( x ):
        cdef Interval y
        y = Interval()
        mpfi_sin( y.x, x.x )
        return y
    def cos( x ):
        cdef Interval y
        y = Interval()
        mpfi_cos( y.x, x.x )
        return y
    def tan( x ):
        cdef Interval y
        y = Interval()
        mpfi_tan( y.x, x.x )
        return y
    def asin( x ):
        cdef Interval y
        y = Interval()
        mpfi_asin( y.x, x.x )
        return y
    def acos( x ):
        cdef Interval y
        y = Interval()
        mpfi_acos( y.x, x.x )
        return y
    def atan( x ):
        cdef Interval y
        y = Interval()
        mpfi_atan( y.x, x.x )
        return y
    def sinh( x ):
        cdef Interval y
        y = Interval()
        mpfi_sinh( y.x, x.x )
        return y
    def cosh( x ):
        cdef Interval y
        y = Interval()
        mpfi_cosh( y.x, x.x )
        return y
    def tanh( x ):
        cdef Interval y
        y = Interval()
        mpfi_tanh( y.x, x.x )
        return y
    def asinh( x ):
        cdef Interval y
        y = Interval()
        mpfi_asinh( y.x, x.x )
        return y
    def acosh( x ):
        cdef Interval y
        y = Interval()
        mpfi_acosh( y.x, x.x )
        return y
    def atanh( x ):
        cdef Interval y
        y = Interval()
        mpfi_atanh( y.x, x.x )
        return y

def from_string(s):
    # unfortunately, we lose precision going through mpf type.
    assert s.startswith('Interval(')
    s = s[len('Interval('):-1]
    l, u = s.split(',')
    l = mpf(l)
    u = mpf(u)
    return Interval(l,u)
    
