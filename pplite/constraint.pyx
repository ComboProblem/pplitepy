# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint


from gmpy2 cimport GMPy_MPZ_From_mpz, import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from cython.operator cimport dereference as deref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from .linear_algebra import Variable, Affine_Expression, Linear_Expression
# from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer
# from .integer_conversions import FLINT_Integer_Conversion_Check
# Using imported integer conversions breaks the code.


import_gmpy2()

cdef class Constraint(object):
    r"""
    Wrapper for PPLite's ``Constraint`` class.

    An object of the class ``Constraint`` is either:

    * an equality :math:`\sum_{i=0}^{n-1} a_i x_i + b = 0`

    * a non-strict inequality :math:`\sum_{i=0}^{n-1} a_i x_i + b \geq 0`

    * a strict inequality :math:`\sum_{i=0}^{n-1} a_i x_i + b > 0`

    where :math:`n` is the dimension of the space, :math:`a_i` is the integer
    coefficient of variable :math:`x_i`, and :math:`b_i` is the integer
    inhomogeneous term.

    INPUT/OUTPUT:

    You construct constraints by writing inequalities in :class:`Variable`,
    :class:`Linear_Expression`, or :class:`Affine_Expression`. Do not attempt to manually construct
    constraints.

    Examples:

    >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> 5*x-2*y+2 > x+y-1
    4*x0-3*x1+3>0
    >>> 5*x-2*y+1 >= x+y
    4*x0-3*x1+1>=0
    >>> 5*x-2*y == x+y
    4*x0-3*x1==0
    >>> 5*x-2*y <= x+y
    -4*x0+3*x1>=0
    >>> 5*x-2*y <  x+y
    -4*x0+3*x1>0
    >>> x+y > 0
    x0+x1>0
    >>> x <= 0
    -x0>=0
    >>> constraint_list = [x+y >= 2, x>0, y<1] 
    >>> constraint_list
    [x0+x1-2>=0, x0>0, -x1+1>0]

    Special care is needed if the left hand side is a constant:

    >>> 0 == 1    # watch out!
    False
    >>> Affine_Expression(0) == Affine_Expression(1)
    -1==0
    >>> Affine_Expression(0) == 1
    -1==0
    """
    def __init__(self, arg=None):
        if arg is None:
            self.thisptr = new Con()
        elif isinstance(arg, Constraint):
            self.thisptr = new Con((<Constraint> arg).thisptr[0])
        else:
            raise TypeError("invalid argument for Constraint")

    def __cinit__(self):
        """
        The Cython constructor.

        See :class:`Constraint` for documentation.

        Tests:

            # >>> from pplite import Variable
            # >>> x = Variable(0)
            # >>> x>0   # indirect doctest
            # x0>0
        """
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        raise TypeError('Constraint unhashable')

    def __repr__(self):
        s = repr(self.linear_form() + self.inhomogeneous_term())
        if self.thisptr.is_equality():
            s += '==0'
        elif self.thisptr.is_strict_inequality():
            s += '>0'
        elif self.thisptr.is_nonstrict_inequality():
            s += '>=0'
        else:
            raise RuntimeError
        return s

    def type(self):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+1== x+y
        >>> c.type()
        'equality'
        >>> c = 3*x+2*y-2<=-x
        >>> c.type()
        'nonstrict_inequality'
        >>> c = x-y > 1
        >>> c.type()
        'strict_inequality'
        """
        cdef ConType t
        t = self.thisptr.type()
        if t == ConType.EQUALITY: 
            return 'equality'
        if t == ConType.NONSTRICT_INEQUALITY:
            return 'nonstrict_inequality'
        if t == ConType.STRICT_INEQUALITY: # ???? Check in on this.
            return 'strict_inequality'
        else:
            raise RuntimeError

    def coefficient(self, v):
        r"""
        INPUT: ``:class:Variable``

        OUTPUT: The coefficient of the ``Variable``.
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+1>= x+y
        >>> c.coefficient(x)
        mpz(4)
        >>> c.coefficient(0)
        mpz(4)
        >>> c.coefficient(1)
        mpz(-3)
        """
        cdef Var* vv
        if isinstance(v, Variable):
            vv = (<Variable> v).thisptr
        else:
            var = Variable(v)
            vv = (<Variable> var).thisptr
        # else:
        #     raise ValueError("Input is not a variable or an integer convertible to a FLINT_Integer.")
        cdef FLINT_Integer coeff
        coeff = self.thisptr.coeff(vv[0])
        return FLINT_Integer_to_Python(coeff)

    def coefficients(self):
        r"""
        INPUT: None
        OUPUT: A tuple of coefficients 
        """
        

    def linear_form(self):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+1>= x+y
        >>> c.linear_form()
        4*x0-3*x1
        """
        cdef Linear_Expr e
        e = self.thisptr.linear_expr()
        linear_form = Linear_Expression()
        linear_form.thisptr[0] = e
        return linear_form

    def inhomogeneous_term(self):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+1>= x+y
        >>> c.inhomogeneous_term()
        mpz(1)
        """
        cdef FLINT_Integer inhom
        inhom = self.thisptr.inhomo_term()
        return FLINT_Integer_to_Python(inhom)

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT: Integer
        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x>=0).space_dimension()
        1
        >>> (y==1).space_dimension()
        2
        """
        return self.thisptr.space_dim()

    def is_equality(self):
        """
        Returns true if the constraint is an equality. 

        OUTPUT: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+2 == x+y-1
        >>> c.is_equality()
        True
        >>> c_2 = 3*x+2*y >= x
        >>> c_2.is_equality()
        False
        >>> c_3 = 3*x-3*y > 1
        >>> c_3.is_equality()
        False
        """
        return self.thisptr.is_equality()

    def is_inequality(self):
        """
        Returns true if the constraint is an inequality (<=, >=, <, >). 

        OUTPUT: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+2 == x+y-1
        >>> c.is_inequality()
        False
        >>> c_2 = 3*x+2*y >= x
        >>> c_2.is_inequality()
        True
        >>> c_3 = 3*x-3*y > 1
        >>> c_3.is_inequality()
        True
        """
        return self.thisptr.is_inequality()

    def is_nonstrict_inequality(self):
        """
        Returns true if the constraint is an nonstrict inequality (<=, >=). 

        OUTPUT: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+2 <= x+y-1
        >>> c.is_nonstrict_inequality()
        True
        >>> c_2 = 3*x+2*y == x
        >>> c_2.is_nonstrict_inequality()
        False
        >>> c_3 = 3*x-3*y < 1
        >>> c_3.is_nonstrict_inequality()
        False
        """
        return self.thisptr.is_nonstrict_inequality()

    def is_strict_inequality(self):
        """
        Returns true if the constraint is a strict inequality (<, >). 

        OUTPUT: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> c = 5*x-2*y+2 > x+y-1
        >>> c.is_strict_inequality()
        True
        >>> c_2 = 3*x+2*y >= x
        >>> c_2.is_strict_inequality()
        False
        >>> c_3 = 3*x-3*y == 1
        >>> c_3.is_strict_inequality()
        False        
        """
        return self.thisptr.is_strict_inequality()

    def zero_dim_false(self):
        """
        TESTS:
        """
        c = Constraint()
        c.thisptr[0] = self.thisptr.zero_dim_false()
        return c

    def zero_dim_positivity(self):
        """
        TESTS:
        """
        c = Constraint()
        c.thisptr[0] = self.thisptr.zero_dim_positivity()
        return c

    def is_tautological(self):
        """
        Tests:
        """
        return self.thisptr.is_tautological()

    def is_inconsistent(self):
        return self.thisptr.is_inconsistent()

    def is_equal_to(self, other):
        cdef Con y
        if isinstance(other, Constraint):
            y = (<Constraint> other).thisptr[0]
            return self.thisptr.is_equal_to(y)
        raise ValueError("This method only works to compare two constraints.")

    def check_inv(self):
        return self.thisptr.check_inv()

    def m_swap(self, other):
        cdef Con y
        if isinstance(other, Constraint):
            y = (<Constraint> other).thisptr[0]
            self.thisptr.m_swap(y)
        raise ValueError("Input must be a Constraint.")

    def set_type(self, t):
        cdef ConType tt
        if t == 'equality': 
            tt = ConType.EQUALITY
            self.thisptr.set_type(tt)
        if t == 'nonstrict_inequality':
            tt = ConType.NONSTRICT_INEQUALITY
            self.thisptr.set_type(tt)
        if t == 'strict_inequality':
            tt = ConType.STRICT_INEQUALITY
            self.thisptr.set_type(tt)
        else:
            raise RuntimeError      

    def is_line_or_equality(self):
        return self.thisptr.is_line_or_equality()

    def set_is_line_or_equality(self):
        self.thisptr.set_is_line_or_equality()

    def linear_combine(self, y, dim):
        cdef Con yy
        cdef dim_type ddim = dim
        if isinstance(y, Constraint) and FLINT_Integer_Conversion_Check(dim):
            yy =  (<Constraint> y).thisptr[0]
            # ddim = Python_int_to_FLINT_Integer(dim)
            self.thisptr.linear_combine(yy, ddim)
        raise ValueError("Types are not correct for this function.")

    def sign_normalize(self):
        self.thisptr.sign_normalize()

    def strong_normalize(self):
        self.thisptr.strong_normalize()

    def check_strong_normalized(self):
        return self.thisptr.check_strong_normalized()


cdef _make_Constraint_from_richcmp(lhs_, rhs_, op):
    cdef Affine_Expr lhs_a
    cdef Affine_Expr rhs_a
    cdef Linear_Expr lhs_e
    cdef Linear_Expr rhs_e
    cdef Var* lhs_v
    cdef Var* rhs_v
    cdef FLINT_Integer lhs_int
    cdef FLINT_Integer rhs_int
    cdef Con constraint
    if isinstance(lhs_, Affine_Expression) and isinstance(rhs_, Affine_Expression):    
        lhs_a = (<Affine_Expression> lhs_).thisptr[0]   
        rhs_a = (<Affine_Expression> rhs_).thisptr[0]        
        if op == Py_LT:
            return _wrap_Constraint(lhs_a < rhs_a)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_a <= rhs_a)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_a == rhs_a)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_a > rhs_a)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_a >= rhs_a)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)   
    if isinstance(lhs_, Linear_Expression) and isinstance(rhs_, Linear_Expression):
        lhs_e = (<Linear_Expression> lhs_).thisptr[0]
        rhs_e = (<Linear_Expression> rhs_).thisptr[0]        
        if op == Py_LT:
            return _wrap_Constraint(lhs_e < rhs_e)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_e <= rhs_e)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_e == rhs_e)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_e > rhs_e)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_e >= rhs_e)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Variable) and isinstance(rhs_, Affine_Expression):    
        lhs_v = (<Variable> lhs_).thisptr
        rhs_a = (<Affine_Expression> rhs_).thisptr[0]        
        if op == Py_LT:
            return _wrap_Constraint(lhs_v[0] < rhs_a)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_v[0] <= rhs_a)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_v[0] == rhs_a)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_v[0] > rhs_a)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_v[0] >= rhs_a)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Affine_Expression) and isinstance(rhs_, Variable):    
        lhs_a = (<Affine_Expression> lhs_).thisptr[0]
        rhs_v = (<Variable> rhs_).thisptr        
        if op == Py_LT:
            return _wrap_Constraint(lhs_a < rhs_v[0])
        elif op == Py_LE:
            return _wrap_Constraint(lhs_a <= rhs_v[0])
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_a == rhs_v[0])
        elif op == Py_GT:
            return _wrap_Constraint(lhs_a > rhs_v[0])
        elif op == Py_GE:
            return _wrap_Constraint(lhs_a >= rhs_v[0])
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Variable) and isinstance(rhs_, Variable):    
        lhs_v = (<Variable> lhs_).thisptr
        rhs_v = (<Variable> rhs_).thisptr        
        if op == Py_LT:
            return _wrap_Constraint(lhs_v[0] < rhs_v[0])
        elif op == Py_LE:
            return _wrap_Constraint(lhs_v[0] <= rhs_v[0])
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_v[0] == rhs_v[0])
        elif op == Py_GT:
            return _wrap_Constraint(lhs_v[0] > rhs_v[0])
        elif op == Py_GE:
            return _wrap_Constraint(lhs_v[0] >= rhs_v[0])
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Affine_Expression) and isinstance(rhs_, Linear_Expression):
        lhs_a = (<Affine_Expression> lhs_).thisptr[0]
        rhs_a = (<Affine_Expression> Affine_Expression(rhs_, 0)).thisptr[0]        
        if op == Py_LT:
            return _wrap_Constraint(lhs_a < rhs_a)
        elif op == Py_LE:           
            return _wrap_Constraint(lhs_a <= rhs_a)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_a == rhs_a)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_a > rhs_a)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_a >= rhs_a)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(rhs_, Affine_Expression) and isinstance(lhs_, Linear_Expression):
        rhs_a = (<Affine_Expression> lhs_).thisptr[0]
        lhs_a = (<Affine_Expression> Affine_Expression(rhs_, 0)).thisptr[0]        
        if op == Py_LT:
            return _wrap_Constraint(lhs_a < rhs_a)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_a <= rhs_a)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_a == rhs_a)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_a > rhs_a)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_a >= rhs_a)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)   
    if isinstance(lhs_, Affine_Expression) and FLINT_Integer_Conversion_Check(rhs_): 
        lhs_a = (<Affine_Expression> lhs_).thisptr[0]
        rhs_int = Python_int_to_FLINT_Integer(rhs_) 
        if op == Py_LT:
            return _wrap_Constraint(lhs_a < rhs_int)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_a <= rhs_int)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_a == rhs_int)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_a > rhs_int)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_a >= rhs_int)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(rhs_, Affine_Expression) and FLINT_Integer_Conversion_Check(lhs_):    
        rhs_a = (<Affine_Expression> lhs_).thisptr[0]
        lhs_int = Python_int_to_FLINT_Integer(rhs_)        
        if op == Py_LT:
            return _wrap_Constraint(lhs_int < rhs_a)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_int <= rhs_a)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_int == rhs_a)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_int > rhs_a)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_int >= rhs_a)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(rhs_, Linear_Expression) and FLINT_Integer_Conversion_Check(lhs_):    
        rhs_e = (<Linear_Expression> rhs_).thisptr[0]
        lhs_int = Python_int_to_FLINT_Integer(lhs_)        
        if op == Py_LT:
            return _wrap_Constraint(lhs_int < rhs_e)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_int <= rhs_e)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_int == rhs_e)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_int > rhs_e)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_int >= rhs_e)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Linear_Expression) and FLINT_Integer_Conversion_Check(rhs_):
        lhs_e = (<Linear_Expression> lhs_).thisptr[0]
        rhs_int = Python_int_to_FLINT_Integer(rhs_)
        if op == Py_LT:
            return _wrap_Constraint(lhs_e < rhs_int)
        elif op == Py_LE:
            return _wrap_Constraint(lhs_e <= rhs_int)
        elif op == Py_EQ:
            return _wrap_Constraint(lhs_e == rhs_int)
        elif op == Py_GT:
            return _wrap_Constraint(lhs_e > rhs_int)
        elif op == Py_GE:
            return _wrap_Constraint(lhs_e >= rhs_int)
        elif op == Py_NE:
            raise NotImplementedError
        else:
            assert(False)
    if isinstance(lhs_, Variable): # variable not explitly declared. promote var to linear expression and compare. 
        lhs_ = Linear_Expression(lhs_)
        return _make_Constraint_from_richcmp(lhs_, rhs_, op)
    if isinstance(rhs_, Variable): # variable not explitly declared. promote var to linear expression and compare. 
        rhs_ = Linear_Expression(rhs_)
        return _make_Constraint_from_richcmp(lhs_, rhs_, op)


cdef _wrap_Constraint(Con constraint):
    """
    Wrap a C++ ``Con`` into a Cython ``Constraint``.
    """
    wrapped_constraint = Constraint()
    wrapped_constraint.thisptr[0] = constraint
    return wrapped_constraint


# Reproduction of these functions here is necessary. Removing this causes things to break.
# Investigation is pending.
cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
    r""" Converts FLINT_Integer to python integer."""
    cdef mpz_t new_int
    mpz_init(new_int)
    fmpz_get_mpz(new_int, integer.impl())
    y = GMPy_MPZ_From_mpz(new_int)
    mpz_clear(new_int)
    return y


cdef FLINT_Integer Python_int_to_FLINT_Integer(integer):
    cdef fmpz_t x
    cdef fmpz y
    if isinstance(integer, (int, str)):
        fmpz_init(x)
        fmpz_set_si(x, integer)
        return FLINT_Integer(x)
    if MPZ_Check(integer): # is this okay?
        y = <fmpz> integer
        return FLINT_Integer(y)
    raise ValueError("Integer Conversion Failed")


def FLINT_Integer_Conversion_Check(possible_integer):
    """
    Checks a python object is convertable to a FLINT_Integer.

    Input: Object

    Output: Bool
    """
    if isinstance(possible_integer, (int, str)):
        return True
    if MPZ_Check(possible_integer):
        return True
    return False