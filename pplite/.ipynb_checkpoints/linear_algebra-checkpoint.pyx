# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from libcpp.vector cimport vector as cppvector

import_gmpy2()

# helper functions for object conversion. 
# It is assumed that the pplite enviroment is set up to use FLINT_integers.
# TODO:  Write a proper conversion module to handle the Integer class in PPLite so this works regardless of setup.

cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
    r""" Converts FLINT_Integer to python integer."""
    
 #    TESTS::
 #        >>> import pplite
 #        >>> cdef fmpz_t x
 #        >>> fmpz_init(x)
 #        >>> fmpz_set_si(x, 7)
 #        >>> cdef FLINT_Integer& w
 #        >>> w = new FLINT_Integer(x)
 #        >>> fmpz_clear(x)
 #        >>> z = FLINT_Integer_to_Python(w)
 #        >>> print(z)
 #        7
 # #   """
    cdef mpz_t new_int
    mpz_init(new_int)
    fmpz_get_mpz(new_int, integer.impl())
    y = GMPy_MPZ_From_mpz(new_int)
    mpz_clear(new_int)
    return y

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer):
    cdef fmpz_t x
    fmpz_init(x)
    if isinstance(integer, (int, str)):
        fmpz_set_si(x, integer)
    return FLINT_Integer(x)

@cython.freelist(128)
cdef class Variable(object):
    r"""
    Wrapper for PPL's ``Variable`` class.

    A dimension of the vector space.    
    
    OTHER STUFF HERE

    INPUT:

    - ``i`` -- integer. The index of the axis.

    OUTPUT:

    A :class:`Variable`

    Examples:

    >>> from pplite import Variable
    >>> x = Variable(123)
    >>> x.id()
    123
    >>> x
    x123

    Note that the "meaning" of an object of the class Variable is completely
    specified by the integer index provided to its constructor: be careful not
    to be mislead by C++ language variable names. For instance, in the following
    example the linear expressions ``e1`` and ``e2`` are equivalent, since the
    two variables ``x`` and ``z`` denote the same Cartesian axis:

    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> z = Variable(0)
    >>> e1 = x + y; e1
    x0+x1
    >>> e2 = y + z; e2
    x0+x1
    >>> e1 - e2
    0
    """
    def __cinit__(self, dim_type i):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        Tests:

        >>> from pplite import Variable
        >>> Variable(123)   # indirect doctest
        x123
        """
        self.thisptr = new Var(i)

        
    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr
        
    def __hash__(self):
        r"""
        Tests:

        >>> import pplite
        >>> hash(pplite.Variable(12))
        Traceback (most recent call last):
        ...
        TypeError: Variable unhashable
        """
        raise TypeError('Variable unhashable')

    def id(self):
        """
        Return the index of the Cartesian axis associated to the variable.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(123)
        >>> x.id()
        123
        """
        return self.thisptr.id()
        
    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer. The returned value is ``self.id()+1``.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> x.space_dimension()
        1
        """
        return self.thisptr.space_dim()
    
    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> x.__repr__()
        'x0'
        """
        return 'x{0}'.format(self.id())
    
    def __add__(self, other):
        r"""
        Return the sum ``self`` + ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` + ``other``.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0);
        >>> y = Variable(1)
        >>> x + y
        x0+x1

        """
        return Linear_Expression(self) + Linear_Expression(other)
        # >>> 15 + y
        # x1+15

        # >>> from gmpy2 import mpz
        # >>> x + mpz(3)
        # x0+3
        # >>> mpz(-5) + y
        # x1-5

        # >>> x + 1.5
        # Traceback (most recent call last):
        # ...
        # TypeError: pplite coefficients must be integral
        # >>> 1.5 + x
        # Traceback (most recent call last):
        # ...
        # TypeError: pplite coefficients must be integral
    def __radd__(self, other):
        return Linear_Expression(self) + Linear_Expression(other)

    def __sub__(self, other):
        r"""
        Return the difference ``self`` - ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` - ``other``.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x - y
        x0-x1
        """
        return Linear_Expression(self) - Linear_Expression(other)

    def __rsub__(self, other):
        return Linear_Expression(other) - Linear_Expression(self)

    def __mul__(self, other):
        r"""
        Return the product ``self`` * ``other``.

        INPUT:

        - ``self``, ``other`` -- One must be an integer, the other a
          :class:`Variable`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` * ``other``.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x * 15
        15*x0
        >>> 15 * y
        15*x1
        """
        #         >>> 1.5 * x
        # Traceback (most recent call last):
        # ...
        # TypeError: pplite coefficients must be integral
        # >>> x * 1.5
        # Traceback (most recent call last):
        # ...
        # TypeError: pplite coefficients must be integral
        if isinstance(self, Variable):
            return Linear_Expression(self) * other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            return Linear_Expression(other) * self

    def __rmul__(self, other):
        return Linear_Expression(self) * other

    def __pos__(self):
        r"""
        Return ``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``+self``

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0); x
        x0
        >>> +x
        x0
        """
        return Linear_Expression(self)

    def __neg__(self):
        r"""
        Return -``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``-self``

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0); x
        x0
        >>> -x
        -x0
        """
        return Linear_Expression(self)*(-1)

    # def __richcmp__(self, other, op):
    #     """
    #     Construct :class:`Constraint` from equalities or inequalities.

    #     INPUT:

    #     - ``self``, ``other`` -- anything convertible to a
    #       :class:`Linear_Expression`

    #     - ``op`` -- the operation.

    #     Examples:

    #     >>> from pplite import Variable
    #     >>> x = Variable(0)
    #     >>> y = Variable(1)
    #     >>> x <  y
    #     -x0+x1>0
    #     >>> x <= 0
    #     -x0>=0
    #     >>> x == y-y
    #     x0==0
    #     >>> x >= -2
    #     x0+2>=0
    #     >>> x >  0
    #     x0>0
    #     >>> 0 == 1    # watch out!
    #     False
    #     >>> 0*x == 1
    #     -1==0
    #     """
    #     return _make_Constraint_from_richcmp(self, other, op)

####################################################
### Linear_Expression ##############################
####################################################
cdef class Linear_Expression(object):
    r"""
    Wrapper for PPL's ``PPL_Linear_Expression`` class.
    """
    def __init__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        """
        cdef size_t i
        if len(args) == 2:
            a = args[0]
            b = args[1]
            self.thisptr = new Linear_Expr()
            if isinstance(a, dict):
                if a:
                    self.thisptr.set_space_dim(1 + max(a))
                    for i, coeff in a.items():
                        self.thisptr.impl()[Variable(i).id()] = Python_int_to_FLINT_Integer(coeff)
            else:
                self.thisptr.set_space_dim(len(a))
                for i, coeff in enumerate(a):
                    self.thisptr.impl()[Variable(i).id()] =  Python_int_to_FLINT_Integer(coeff)
            return 
            # self.thisptr.set_inhomogeneous_term(PPL_Coefficient_from_pyobject(b))
            # return
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, Variable):
                v = <Variable>arg
                self.thisptr = new Linear_Expr(v.thisptr[0])
                return
            if isinstance(arg, Linear_Expression):
                e = <Linear_Expression>arg
                self.thisptr = new Linear_Expr(e.thisptr[0])
                return
            #self.thisptr. 
            #raise ValueError("Initalizing with one argument requires either a linear expression or a variable to be passed in.")
        elif len(args) == 0:
            self.thisptr = new Linear_Expr()
            return
        else:
            raise ValueError("Cannot initialize with more than 2 arguments.")

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import pplite
        >>> hash(pplite.Linear_Expression(10))
        Traceback (most recent call last):
        ...
        TypeError: Linear_Expression unhashable
        """
        raise TypeError('Linear_Expression unhashable')

    def space_dimension(self):
        """
        Return the dimension of the vector space necessary for the
        linear expression.

        OUTPUT:

        Integer.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x+y).space_dimension()
        2
        >>> (x+y).space_dimension()
        2
        >>> (y).space_dimension()
        2
        >>> (x).space_dimension()
        1
        >>> (y-y).space_dimension()
        2
        """
        return self.thisptr.space_dim()
    
    def set_space_dimension(self, dim_type dim):
        self.thisptr.set_space_dim(dim)
        
    def coefficient(self, v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An (Python) Integer. 

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> e = 3*x
        >>> e.coefficient(x)
        mpz(3)
        """
        #      >>> e.coefficient(Variable(1))
        # mpz(0)   
        cdef Variable vv # rewrite this method to read coeffs correctly
        
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        return FLINT_Integer_to_Python(self.thisptr.impl()[vv.id()])
    
    def set_coefficient(self, i, v):
        """
        Set the ``i``-th coefficient to ``v``.

        INPUT:

        - ``i`` - variable or variable index

        - ``v`` - integer

        Examples:

        >>> from pplite import Variable
        >>> L = Variable(0) + 3 * Variable(1)
        >>> L.set_coefficient(1, -5)
        >>> L.set_coefficient(7, 3)
        >>> L
        x0-5*x1
        """
        cdef Variable ii
        if type(i) is Variable:
            ii = <Variable> i
        else:
            ii = Variable(i)
        cdef FLINT_Integer vv 
        vv = Python_int_to_FLINT_Integer(v)
        self.thisptr.impl()[ii.id()] = vv
    
    def __repr__(self):
        r"""
        Return a string representation of the linear expression.

        OUTPUT:

        A string.

        Examples:

        >>> from pplite import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x
        x0
        >>> x-x
        0
        >>> 2*x
        2*x0
        """
        s = ''
        first = True
        for i in range(self.space_dimension()):
            x = Variable(i)
            coeff = self.coefficient(x)
            if coeff == 0:
                continue
            if first and coeff == 1:
                s += '%r' % x
                first = False
            elif first and coeff == -1:
                s += '-%r' % x
                first = False
            elif first and coeff != 1:
                s += '%d*%r' % (coeff, x)
                first = False
            elif coeff == 1:
                s += '+%r' % x
            elif coeff == -1:
                s += '-%r' % x
            else:
                s += '%+d*%r' % (coeff, x)
        if first:
            s = '0'
        return s

    def swap_space_dimensions(self, v1, v2):
        r"""
        Swaps the coefficients of ``v1`` and ``v2``.

        INPUT:

        - ``v1``, ``v2`` - variables or indices of variables

        Examples:

        >>> from pplite import Variable
        >>> L = Variable(1) - 3 * Variable(3)
        >>> L.swap_space_dimensions(Variable(1), Variable(3))
        >>> L
        -3*x1+x3

        >>> L = Variable(1) - 3 * Variable(3)
        >>> L.swap_space_dimensions(1, 3)
        >>> L
        -3*x1+x3
        """
        cdef dim_type var_1, var_2
        if type(v1) is Variable:
            var_1 = v1.id()
        else:
            vv1 = new Var(v1)
            var_1 = vv1.id()
        if type(v2) is Variable:
            var_2 = v2.id()
        else:
            vv2 = new Var(v2)
            var_2 = vv2.id()
        self.thisptr.swap_space_dims(var_1, var_2)

    # def shift_space_dimensions(self, v, PPL_dimension_type n):
    #     r"""
    #     Shift by ``n`` the coefficients of variables starting from the
    #     coefficient of ``v``.

    #     This increases the space dimension by ``n``.

    #     Examples:

    #     >>> import ppl
    #     >>> L = pplite.Variable(0) + 13 * pplite.Variable(2) + 5 * pplite.Variable(7)
    #     >>> L
    #     x0+13*x2+5*x7
    #     >>> L.shift_space_dimensions(pplite.Variable(2), 2)
    #     >>> L
    #     x0+13*x4+5*x9
    #     >>> L.shift_space_dimensions(pplite.Variable(7), 3)
    #     >>> L
    #     x0+13*x4+5*x12
    #     """
    #     cdef Variable vv
    #     if type(v) is Variable:
    #         vv = <Variable> v
    #     else:
    #         vv = Variable(v)
    #     self.thisptr.shift_space_dimensions(vv.thisptr[0], n)

    # def remove_space_dimensions(self, Variables_Set V):
    #     r"""
    #     Removes the dimension specified by the set of variables ``V``.

    #     See :class:`Variables_Set` to construct set of variables.

    #     Examples:

    #     >>> import ppl
    #     >>> L = sum(i * pplite.Variable(i) for i in range(10))
    #     >>> L
    #     x1+2*x2+3*x3+4*x4+5*x5+6*x6+7*x7+8*x8+9*x9
    #     >>> L.remove_space_dimensions(pplite.Variables_Set(3,5))
    #     >>> L
    #     x1+2*x2+6*x3+7*x4+8*x5+9*x6
    #     """
    #     self.thisptr.remove_space_dimensions(V.thisptr[0])

    # def all_homogeneous_terms_are_zero(self):
    #     """
    #     Test if ``self`` is a constant linear expression.

    #     OUTPUT:

    #     Boolean.

    #     Examples:

    #     >>> from pplite import Variable, Linear_Expression
    #     >>> Linear_Expression(10).all_homogeneous_terms_are_zero()
    #     True
    #     """
    #     return self.thisptr.all_homogeneous_terms_are_zero()

    def is_equal_to(self, Linear_Expression other):
        """
        Test equality with another linear expression.

        OUTPUT: boolean

        Examples:

        >>> from pplite import Variable
        >>> L1 = Variable(0) + 2 * Variable(3)
        >>> L2 = Variable(0) + 2 * Variable(3)
        >>> L3 = Variable(0) - Variable(2)
        >>> L1.is_equal_to(L2)
        True
        >>> L1.is_equal_to(L3)
        False
        """
        return self.thisptr.is_equal_to(other.thisptr[0])

    # def ascii_dump(self):
    #     r"""
    #     Write an ASCII dump to stderr.

    #     Examples:

    #     >>> cmd  = 'from pplite import Linear_Expression, Variable\n'
    #     >>> cmd += 'x = Variable(0)\n'
    #     >>> cmd += 'y = Variable(1)\n'
    #     >>> cmd += 'e = 3*x+2*y+1\n'
    #     >>> cmd += 'e.ascii_dump()\n'
    #     >>> from subprocess import Popen, PIPE
    #     >>> import sys
    #     >>> proc = Popen([sys.executable, '-c', cmd], stdout=PIPE, stderr=PIPE)
    #     >>> out, err = proc.communicate()
    #     >>> len(out) == 0
    #     True
    #     >>> len(err) > 0
    #     True
    #     """
    #     self.thisptr.ascii_dump()
        
    def __add__(self, other):
        r"""
        Add ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Linear_Expression`

        Examples:

        >>> from pplite import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x + y + y + y
        x0+3*x1
        """
        # >>> from gmpy2 import mpz
        # >>> mpz(3) + x + mpz(5) + y + mpz(7)
        # x0+x1+15
        cdef Linear_Expr* lhs
        cdef Linear_Expr* rhs
        
        if isinstance(self, Linear_Expression):
            lhs = (<Linear_Expression> self).thisptr
        else:
            lhs_expr = Linear_Expression(self)
            lhs = (<Linear_Expression> lhs_expr).thisptr
        if isinstance(other, Linear_Expression):
            rhs = (<Linear_Expression> other).thisptr
        else:
            rhs_expr = Linear_Expression(other)
            rhs = (<Linear_Expression> rhs_expr).thisptr
        cdef Linear_Expr result
        result = lhs[0] + rhs[0]
        result_expr = Linear_Expression()
        # result_expr
        result_expr.thisptr[0] = result #could be copying or moving?
        return result_expr

    def __radd__(self, other):
        cdef Linear_Expr* lhs
        cdef Linear_Expr* rhs

        lhs = (<Linear_Expression> self).thisptr
        if isinstance(other, Linear_Expression):
            rhs = (<Linear_Expression> other).thisptr
        else:
            rhs_expr = Linear_Expression(other)
            rhs = (<Linear_Expression> rhs_expr).thisptr
        cdef Linear_Expr result
        result = lhs[0] + rhs[0]
        result_expr = Linear_Expression()
        # result_expr
        result_expr.thisptr[0] = result
        return result_expr

    def __sub___(self, other):
        cdef Linear_Expr* lhs
        cdef Linear_Expr* rhs
        
        if isinstance(self, Linear_Expression):
            lhs = (<Linear_Expression> self).thisptr
        else:
            lhs_expr = Linear_Expression(self)
            lhs = (<Linear_Expression> lhs_expr).thisptr
        if isinstance(other, Linear_Expression):
            rhs = (<Linear_Expression> other).thisptr
        else:
            rhs_expr = Linear_Expression(other)
            rhs = (<Linear_Expression> rhs_expr).thisptr
        cdef Linear_Expr result
        result = lhs[0] - rhs[0]
        result_expr = Linear_Expression()
        # result_expr
        result_expr.thisptr[0] = result #could be copying or moving?
        return result_expr

    def __rsub__(self, other):
        cdef Linear_Expr* lhs
        cdef Linear_Expr* rhs

        lhs = (<Linear_Expression> self).thisptr
        if isinstance(other, Linear_Expression):
            rhs = (<Linear_Expression> other).thisptr
        else:
            rhs_expr = Linear_Expression(other)
            rhs = (<Linear_Expression> rhs_expr).thisptr
        cdef Linear_Expr result
        result = rhs[0] - lhs[0]
        result_expr = Linear_Expression()
        # result_expr
        result_expr.thisptr[0] = result
        return result_expr



    def __mul__(self, other):
        r"""
        Multiply ``self`` with ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The product as a :class:`Linear_Expression`

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> 8 * (x)
        8*x0
        >>> y * 8
        8*x1
        """
        #         >>> 2**128 * x
        # 340282366920938463463374607431768211456*x0
        #         >>> from gmpy2 import mpz
        # >>> mpz(3) * x * mpz(5)
        # 15*x0
        cdef Linear_Expr* e
        
        if isinstance(self, Linear_Expression):
            e = (<Linear_Expression> self).thisptr
            c = other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            e = (<Linear_Expression> other).thisptr
            c = self

        cdef FLINT_Integer cc = Python_int_to_FLINT_Integer(c)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = e[0] * cc
        return result


def test_current_obj():
    x = Variable(0)
    e = Linear_Expression(x)
    x_2 = Variable(2)
    # print(e.coefficient(x))
    # print(e.coefficient(x_2))
    # cdef FLINT_Integer w
    # w = Python_int_to_FLINT_Integer(2)
    # print(FLINT_Integer_to_Python(w))
    e.set_coefficient(x_2, 4)
    # print(e.coefficient(x_2))
    e_2 = Linear_Expression([1, 2, 3, 4], 5)
    print(e_2, "e_2")
    print(e, "e")
    print(e + e_2, "e_2 + e")
    print(Variable(1) + e_2, "x1+e_2")
    print(e_2*4, "e_2*4")
    # print(e_2)
    # e_2.swap_space_dimensions(Variable(2),Variable(3))
    # print(e_2)
    # print(Variable(2), Variable(3))