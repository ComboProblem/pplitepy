# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from libcpp.vector cimport vector as cppvector
from .constraint cimport _make_Constraint_from_richcmp
# from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer

# TODO: Investigate why everything breaks when importing integer conversion as opposed to local definitions. 

import_gmpy2()

cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
    r""" Converts PPLite::FLINT_Integer to python integer."""
    cdef mpz_t new_int
    mpz_init(new_int)
    fmpz_get_mpz(new_int, integer.impl())
    y = GMPy_MPZ_From_mpz(new_int)
    mpz_clear(new_int)
    return y

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer):
    r"""Converts python sting or int to a PPLite::FLINT_Integer."""
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

################
### Variable ###
################

@cython.freelist(256)
cdef class Variable(object):
    r"""
    Wrapper for PPLites's ``Var`` class.

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


    def swap(self, Variable w):
        """
        Swaps representation of variables.
        """
        swap(self.thisptr[0], w.thisptr[0])
        return self
    
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
          ``Affine_Expression``: An integer, a :class:`Variable`,
          :class:`Linear_Expression`, or :class`Affine_Expression`.

        OUTPUT:

        A :class:`Linear_Expression`.

        Examples:

        >>> from pplite import Variable, Linear_Expression, Affine_Expression
        >>> x = Variable(0);
        >>> y = Variable(1)
        >>> x + y
        x0+x1
        >>> isinstance(x+y, Linear_Expression)
        True
        >>> isinstance(x+y, Affine_Expression)
        False
        """
        if not isinstance(self, Variable):
            raise NotImplemented 
        if not isinstance(other, Variable): # promote variable to linear expression. 
            return Linear_Expression(self) +  other
        # only use this method if both classes are an instance of Variable
        other_var = <Variable> other
        result = Linear_Expression()
        result.thisptr = new Linear_Expr(self.thisptr[0] + other_var.thisptr[0])
        return result

    def __radd__(self, other):
        """
        Return the sum ``other`` + ``self``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Affine_Expression``: An integer, a :class:`Variable`,
          :class:`Linear_Expression`, or :class`Affine_Expression`.

        OUTPUT:

        A :class:`Linear_Expression`. 
        """
        if not isinstance(self, Variable):
            raise NotImplemented
        if not isinstance(other, Variable):
            return other + Linear_Expression(self)
        other_var = <Variable> other
        result = Linear_Expression()
        result.thisptr = new Linear_Expr(self.thisptr[0] + other_var.thisptr[0])
        return result

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
        if not isinstance(self, Variable):
            return NotImplemented
        if not isinstance(other, Variable):
            return Linear_Expression(self) - other
        other_var = <Variable> other
        result = Linear_Expression()
        result.thisptr = new Linear_Expr(self.thisptr[0] - other_var.thisptr[0])
        return result


    def __rsub__(self, other):
        """
        Return the difference ``other`` - ``self``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` - ``other``.
        """
        if not isinstance(self, Variable): 
            return NotImplemented
        if not isinstance(other, Variable):
            return other - Linear_Expression(self)
        other_var = <Variable> other
        result = Linear_Expression()
        result.thisptr = new Linear_Expr(other_var.thisptr[0] - self.thisptr[0])
        return result

    def __mul__(self, other):
        r"""
        Return the product ``self`` * ``other``.

        INPUT:

        - ``self``, ``other`` -- One must be an integer, the other a
          :class:`Variable`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` * ``other``.

        Examples:

        >>> from pplite import Variable, Linear_Expression
        >>> x = Variable(0); y = Variable(1)
        >>> x * 15
        15*x0
        >>> e = 15 * y; e
        15*x1
        >>> isinstance(e, Linear_Expression)
        True
        >>> x * 1.5
        Traceback (most recent call last):
        ...
        ValueError: Integer Conversion Failed
        """
        if isinstance(self, Variable): # pplite doesn't explictly multiply Var and ints, convert to linear expr
            return Linear_Expression(self) * other
        return NotImplemented

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

    def __richcmp__(self, other, op):
        """
        Construct :class:`Constraint` from equalities or inequalities.

        INPUT:

        - ``self``, ``other`` -- anything convertible to a
          :class:`Linear_Expression`

        - ``op`` -- the operation.

        Examples:
        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x <  y
        -x0+x1>0
        >>> x <= 0
        -x0>=0
        >>> x == y-y
        x0==0
        >>> x >= -2
        x0+2>=0
        >>> x >  0
        x0>0
        >>> 0 == 1  # watch out!
        False
        >>> 0*x == 1
        -1==0
        """
        return _make_Constraint_from_richcmp(self, other, op)

#########################
### Linear_Expression ###
#########################

cdef class Linear_Expression(object):
    r"""
    Wrapper for PPLite's ``Linear_Expr`` class.

    This class might more aptly be described as linear form rather than a linear expression. 
    For translation purposes, the class is named linear Expression to align with 
    the original pplite code and ppl.

    INPUT:

    The constructor accepts zero, one, or two arguments.

    If there are two arguments ``Linear_Expression(a,b)``, they are
    interpreted as

    - ``a`` -- a :class:`Linear_Expression`.

    - ``b`` -- an positive integer. The space dimension of a linear form.

    A single argument ``Linear_Expression(arg)`` is interpreted as

    - ``arg`` -- something that determines a linear
      expression. Possibilities are:

      * a :class:`Variable`: The linear expression given by that
        variable.

      * a :class:`Linear_Expression`: The copy constructor.

      * an integer: Constructs the 0 linear expression for space dimension of the integer.

    No argument is the default constructor and returns the zero linear
    expression.

    OUTPUT:

    A :class:`Linear_Expression`

    Examples:

    >>> from pplite import Variable, Linear_Expression

    >>> e = -3*Variable(1) + Variable(7); e
    -3*x1+x7
    >>> isinstance(e, Linear_Expression)
    True
    >>> e.space_dimension()
    8
    >>> e_2 = Linear_Expression(e, 20); e_2
    -3*x1+x7
    >>> e_2.space_dimension()
    20
    >>> Linear_Expression()
    0
    >>> e = Linear_Expression(5); e
    0
    >>> e.space_dimension()
    5
    >>> e = Linear_Expression(Variable(2)); e
    x2
    >>> x = Variable(123)
    >>> y = Variable(321)
    >>> expr = x+y
    >>> expr
    x123+x321
    >>> expr.coefficient(x)
    mpz(1)
    >>> expr.coefficient(Variable(124))
    mpz(0)

    String, rationals and floating point types are accepted as long as they
    represent exact integers:
    """
    # TODO: Finish Tests/claim above.
    def __init__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        """
        cdef dim_type dim 
        if len(args) == 2:
            a = args[0]
            b = args[1]
            dim = b
            if isinstance(a, Linear_Expression):
                e = <Linear_Expression> a
                self.thisptr = new Linear_Expr(e.thisptr[0], b)
            return 
        elif len(args) == 1:
            arg = args[0]
            if isinstance(arg, Variable):
                v = <Variable> arg
                self.thisptr = new Linear_Expr(v.thisptr[0])
                return
            if isinstance(arg, Linear_Expression):
                e = <Linear_Expression> arg
                self.thisptr = new Linear_Expr(e.thisptr[0])
                return
            if isinstance(arg, int):
                dim = arg
                self.thisptr = new Linear_Expr(dim)
                return
            raise ValueError("Initalizing with one argument requires either a linear expression, variable, or integer to be passed in.")
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

        Examples::

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
        """
        Sets the ambient dimension which the linear expression lives in.

        INPUT: 

        - ``dim`` an integer
        """
        self.thisptr.set_space_dim(dim)

    def coefficient(self, v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable` or an integer.

        OUTPUT:

        An (Python) Integer. 

        Examples::

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> e = 3*x
        >>> e.coefficient(x)
        mpz(3)
        >>> e.coefficient(Variable(1))
        mpz(0)
        """
        cdef Var* vv
        cdef FLINT_Integer result
        if isinstance(v, Variable):
            vv = (<Variable> v).thisptr
            result = self.thisptr[0].get(vv[0])
            return FLINT_Integer_to_Python(result)
        cdef dim_type dim
        if isinstance(v, int):
            dim = v
            result = self.thisptr[0].get(dim)
            return FLINT_Integer_to_Python(result)
    
    def set_coefficient(self, i, n):
        """
        Set the ``i``-th coefficient to ``n``.

        INPUT:

        - ``i`` - variable or variable index

        - ``n`` - integer

        Examples::

        >>> from pplite import Variable
        >>> L = Variable(0) + (3 * Variable(1)); L
        x0+3*x1
        >>> L.set_coefficient(1, -5)
        >>> L
        x0-5*x1
        >>> L.set_coefficient(3, 7); L
        x0-5*x1+7*x3
        """
        cdef FLINT_Integer nn
        if isinstance(i, Variable):
            ii = <Variable> i
        else:
            var_i = Variable(i)
            ii = <Variable> var_i
        nn = Python_int_to_FLINT_Integer(n)
        (<Linear_Expression> self).thisptr[0].set(ii.thisptr[0], nn) 

    def __repr__(self):
        r"""
        Return a string representation of the linear expression.

        OUTPUT:

        A string.

        Examples::

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
        >>> L = Variable(1) - 3 * Variable(3); L
        x1-3*x3
        >>> L.swap_space_dimensions(Variable(1), Variable(3))
        >>> L
        -3*x1+x3
        >>> L = Variable(1) - 3 * Variable(3)
        >>> L.swap_space_dimensions(1, 3)
        >>> L
        -3*x1+x3
        """
        cdef dim_type var_1, var_2
        if isinstance(v1, Variable):
            var_1 = v1.id()
        else:
            vv1 = new Var(v1)
            var_1 = vv1.id()
        if isinstance(v2, Variable):
            var_2 = v2.id()
        else:
            vv2 = new Var(v2)
            var_2 = vv2.id()
        self.thisptr.swap_space_dims(var_1, var_2)

    def shift_space_dimensions(self, v, dim_type n):
        r"""
        Shift by ``n`` the coefficients of variables starting from the
        coefficient of ``v``.

        This increases the space dimension by ``n``.

        INPUT:

        - ``v`` a :class:`Variable`

        - ``n`` an integer

        Examples::

        >>> from pplite import Variable
        >>> L = Variable(0) + 13 * Variable(2) + 5 * Variable(7)
        >>> L
        x0+13*x2+5*x7
        >>> L.shift_space_dimensions(Variable(2), 2)
        >>> L
        x0+13*x4+5*x9
        >>> L.shift_space_dimensions(Variable(7), 3)
        >>> L
        x0+13*x4+5*x12
        """
        cdef Variable vv
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        self.thisptr.shift_space_dims(vv.thisptr[0], n)

    # def remove_space_dimensions(self, Variables_Set V):
    #     r"""
    #     Removes the dimension specified by the set of variables ``V``.

    #     See :class:`Variables_Set` to construct set of variables.

    #     Examples:

    #     >>> from pplite import Variable
    #     >>> L = sum(i * Variable(i) for i in range(10))
    #     >>> L
    #     x1+2*x2+3*x3+4*x4+5*x5+6*x6+7*x7+8*x8+9*x9
    #     >>> L.remove_space_dimensions(Variables_Set(3,5))
    #     >>> L
    #     x1+2*x2+6*x3+7*x4+8*x5+9*x6
    #     """
    #     self.thisptr.remove_space_dimensions(V.thisptr[0])

    def all_homogeneous_terms_are_zero(self):
        """
        Test if ``self`` is a constant linear expression.

        OUTPUT:

        Boolean.

        Examples:

        >>> from pplite import Variable, Linear_Expression
        >>> x = Variable(1)
        >>> (x-x).all_homogeneous_terms_are_zero()
        True
        """
        return self.thisptr.is_zero()

    def is_equal_to(self, Linear_Expression other):
        """
        Test equality with another linear expression.

        OUTPUT: boolean

        Examples::

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
        
    def __add__(self, other):
        r"""
        Add ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression` or :class:`Affine_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Linear_Expression` or :class:`Affine_Expression` depending on input.

        Examples::

        >>> from pplite import Affine_Expression, Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x + y + y + y
        x0+3*x1
        >>> e = x + y + x
        >>> isinstance(e, Linear_Expression)
        True
        >>> e2 = 7 + e; e2
        2*x0+x1+7
        >>> isinstance(e2, Linear_Expression)
        False
        >>> isinstance(e2, Affine_Expression)
        True
        >>> e + 7
        2*x0+x1+7
        >>> e + Variable(0)
        3*x0+x1
        >>> Variable(1) + e
        2*x0+2*x1
        >>> e3 = e + e2; e3
        4*x0+2*x1+7
        >>> isinstance(e3, Affine_Expression)
        True
        """
        if not isinstance(self, Linear_Expression):
            return NotImplemented
        # to mimic pplite, we use cases for type conversions.
        # case 1: linear_expr + var -> linear_expr
        if isinstance(other, Variable):
            other_var = <Variable> other
            result = Linear_Expression()
            result.thisptr = new Linear_Expr(other_var.thisptr[0] + self.thisptr[0])  
            return result
        # case 2: linear_expr + linear_expr -> linear_expr
        cdef Linear_Expr* lhs_expr
        cdef Linear_Expr* rhs_expr
        cdef Linear_Expr result_expr
        if isinstance(other, Linear_Expression):
            lhs_expr = (<Linear_Expression> self).thisptr
            rhs_expr = (<Linear_Expression> other).thisptr
            result_expr = lhs_expr[0] + rhs_expr[0]
            result = Linear_Expression()
            result.thisptr[0] = result_expr
            return result
        # case 3: linear_expr + affine_expr -> affine_expr
        if isinstance(other, Affine_Expression):
            return NotImplemented 
        # case 4, linear_expr + integer -> affine expression
        result_aff = Affine_Expression(self, other)
        return result_aff

    def __radd__(self, other):
        return self + other        

    def __sub__(self, other):
        r"""
        Subtract ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression` or :class:`Affine_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Linear_Expression` or :class:`Affine_Expression` depending on input.

        Examples:

        >>> from pplite import Affine_Expression, Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x - y + y + y
        x0+x1
        >>> e = y - x + y - x; e
        -2*x0+2*x1
        >>> isinstance(e, Linear_Expression)
        True
        >>> e - e
        0
        >>> e1 = y - x
        >>> e1 - e
        x0-x1
        >>> e2 = e - 7; e2
        -2*x0+2*x1-7
        >>> e3 = 7 - e; e3
        2*x0-2*x1+7
        >>> isinstance(e3, Linear_Expression)
        False
        >>> isinstance(e3, Affine_Expression)
        True
        """
        # to mimic pplite, we use cases for type conversions.
        # case 1: linear_expr - var -> linear_expr
        if not isinstance(self, Linear_Expression):
            return NotImplemented
        if isinstance(other, Variable):
            other_var = <Variable> other
            result = Linear_Expression()
            self.thisptr[0] -= other_var.thisptr[0]
            result.thisptr = new Linear_Expr(self.thisptr[0])  
            return result
        # case 2: linear_expr + linear_expr -> linear_expr
        cdef Linear_Expr* lhs_expr
        cdef Linear_Expr* rhs_expr
        cdef Linear_Expr result_expr
        if isinstance(other, Linear_Expression):
            lhs_expr = (<Linear_Expression> self).thisptr
            rhs_expr = (<Linear_Expression> other).thisptr
            result_expr = lhs_expr[0] - rhs_expr[0]
            result = Linear_Expression()
            result.thisptr[0] = result_expr
            return result
        # case 3: linear_expr + affine_expr -> affine_expr
        if isinstance(other, Affine_Expression):
            return NotImplemented 
        # case 4, linear_expr + integer -> affine expression
        result_aff = Affine_Expression(self, -other)
        return result_aff

    def __rsub__(self, other):
        if not isinstance(self, Linear_Expression):
            return NotImplemented
        if isinstance(other, Variable):
            other_var = <Variable> other
            result = Linear_Expression()
            result.thisptr = new Linear_Expr(other_var.thisptr[0] - self.thisptr[0])  
            return result
        cdef Linear_Expr* lhs_expr
        cdef Linear_Expr* rhs_expr
        cdef Linear_Expr result_expr
        if isinstance(other, Linear_Expression):
            lhs_expr = (<Linear_Expression> self).thisptr
            rhs_expr = (<Linear_Expression> other).thisptr
            result_expr = rhs_expr[0] - lhs_expr[0]
            result = Linear_Expression()
            result.thisptr[0] = result_expr
            return result
        if isinstance(other, Affine_Expression):
            return NotImplemented
        result_aff = Affine_Expression(-1*self, other)
        return result_aff

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

        Examples::

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

    def __rmul__(self, other):
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

    def __neg__(self):
        neg_assign(self.thisptr[0])
        return self

    def __pos__(self):
        return self

    def __richcmp__(self, other, op):
        """
        Construct :class:`Constraint`s

        Examples::

        >>> from pplite import Variable
        >>> A = Variable(0)
        >>> B = Variable(1)
        >>> A > B
        x0-x1>0
        """
        return _make_Constraint_from_richcmp(self, other, op)


#########################
### Affine_Expression ###
#########################

cdef class Affine_Expression(object):
    r"""
    Wrapper for PPLite's ``Affine_Expr`` class.

    Note, this behaves analgolously to ppl's ``Linear_Expression`` class. 

    Examples:

    The constructor accepts zero, one, or two arguments.

    If there are two arguments ``Affine_Expression(a,b)``, they are
    interpreted as

    - ``a`` -- either a dictionary whose indices are space dimension and
      values are coefficients, an iterable coefficients (e.g. a list or
      tuple), or an :class:`Linear_Expression`.

    - ``b`` -- an integer. The inhomogeneous term.

    A single argument ``Affine_Expression(arg)`` is interpreted as

    - ``arg`` -- something that determines a affine
      expression. Possibilities are:

      * a :class:`Affine_Expression`: The copy constructor.

      * an integer: Constructs the constant affine expression.

    No argument is the default constructor and returns the zero affine
    expression.

    OUTPUT:

    A :class:`Affine_Expression`

    Examples::

    >>> from pplite import Variable, Linear_Expression, Affine_Expression
    >>> e = Variable(2) - 3*Variable(4)
    >>> isinstance(e, Linear_Expression)
    True
    >>> a = Affine_Expression(e, 7); a
    x2-3*x4+7
    >>> Affine_Expression()
    0
    >>> -a
    -x2+3*x4-7
    """
    def __init__(self, *args):
        """
        The Cython constructor.

        See :class:`Affine_Expression` for documentation.
        """
        cdef FLINT_Integer k
        if len(args) == 2:
            expr_arg = args[0]
            int_arg = args[1]   
            k = Python_int_to_FLINT_Integer(int_arg)
            if isinstance(expr_arg, Linear_Expression):
                e = <Linear_Expression> expr_arg
                self.thisptr = new Affine_Expr(e.thisptr[0], k)
                return
            raise ValueError("Initializing with two arguments requires a :class:`Linear_Expression` and an integer to be passed in.")
        if len(args) == 1:
            arg = args[0]   
            if isinstance(arg, int):
                k = Python_int_to_FLINT_Integer(arg)
                self.thisptr = new Affine_Expr(k)
                return
            if isinstance(arg, Affine_Expression):
                a = <Affine_Expression> arg 
                self.thisptr = new Affine_Expr(a.thisptr[0])
                return
            raise ValueError("Initializing with one argument requires either an :class:`Affine_Expression` xor an integer to be passed in.")
        elif len(args) == 0:
            self.thisptr = new Affine_Expr()
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
        Tests::

        >>> import pplite
        >>> hash(pplite.Affine_Expression(10))
        Traceback (most recent call last):
        TypeError: Affine_Expression unhashable
        """
        raise TypeError('Affine_Expression unhashable')

    def space_dimension(self):
        """
        Return the dimension of the vector space necessary for the
        linear expression.

        OUTPUT:

        Integer.

        Examples::

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x+y+1).space_dimension()
        2
        """
        return self.thisptr.space_dim()

    def linear_form(self):
        """
        Returns the linear form of the affine expression.

        OUTPUT:

        :class:`Linear_Expression`

        Examples::

        >>> from pplite import Variable, Affine_Expression
        >>> e = Variable(2) - 3*Variable(4)
        >>> a = Affine_Expression(e, 7); a
        x2-3*x4+7
        >>> a.linear_form()
        x2-3*x4
        """
        cdef Linear_Expr e
        e = (<Affine_Expression> self).thisptr[0].expr
        ee = Linear_Expression()
        ee.thisptr[0] = e
        return ee

    def inhomogeneous_term(self):
        """
        Returns the inhogogenous term of an affine expression.

        OUTPUT:

        Integer.

        Examples::

        >>> from pplite import Variable, Affine_Expression
        >>> e = Variable(2) - 3*Variable(4)
        >>> a = Affine_Expression(e, 7); a
        x2-3*x4+7
        >>> a.inhomogeneous_term()
        mpz(7)
        """
        cdef FLINT_Integer c
        c = self.thisptr.inhomo
        return FLINT_Integer_to_Python(c)

    def coefficient(self, v):   
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An Integer. 

        Examples::

        >>> from pplite import Variable, Affine_Expression
        >>> e = Variable(2) - 3*Variable(4)
        >>> a = Affine_Expression(e, 7)
        >>> a.coefficient(Variable(2))
        mpz(1)
        """
        return self.linear_form().coefficient(v)

    def __repr__(self):
        r"""
        Return a string representation of the linear expression.

        OUTPUT:

        A string.

        Examples::

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
        inhomog = self.inhomogeneous_term()
        if inhomog != 0:
            if first:
                s += '%d' % inhomog
                first = False
            else:
                s += '%+d' % inhomog
        if first:
            s = '0'
        return s

    def all_terms_are_zero(self):
        """
        Test if ``self`` is a constant linear expression.

        OUTPUT:

        boolean
        """
        return self.thisptr.is_zero()

    def is_equal_to(self, Affine_Expression other):
        """
        Test equality with another affine expression.

        OUTPUT: 

        boolean
        """
        if self.inhomogeneous_term() == other.inhomogeneous_term() and self.linear_form().is_equal_to(other.linear_form()):
            return True
        return False

    def m_swap(self, Affine_Expression y):
        self.thisptr[0].m_swap(y.thisptr[0])
        return self

    def normalize(self):
        """
        Normalize the Affine_expression.

        """
        self.thisptr[0].normalize()

    def sign_normalize(self):
        self.thisptr[0].sign_normalize()

    def __add__(self, other):
        r"""
        Add ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Affine_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Affine_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Affine_Expression`

        Examples::

        >>> from pplite import Affine_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> a = x + y + y + y + 1; a
        x0+3*x1+1
        >>> a.inhomogeneous_term()
        mpz(1)
        >>> a + a
        2*x0+6*x1+2
        >>> 12 + a
        x0+3*x1+13
        >>> a + Variable(2)
        x0+3*x1+x2+1
        >>> a.linear_form()
        x0+3*x1
        >>> a.linear_form()+ a.inhomogeneous_term()
        x0+3*x1+1
        """
        cdef Affine_Expr* lhs
        if isinstance(self, Affine_Expression):
            lhs = (<Affine_Expression> self).thisptr
        else:
            lhs_expr = Affine_Expression(self)
            lhs = (<Affine_Expression> lhs_expr).thisptr
        result_expr = Affine_Expression()
        cdef Var* vv
        cdef Affine_Expr* temp
        if isinstance(other, Variable):
            other_var = <Variable> other
            vv = other_var.thisptr
            temp = new Affine_Expr(self.thisptr[0])
            self.thisptr[0] += vv[0]  # There is an oddity where adding affine expression and variable doesn't work.
            # This means to add via PPlite add methods, we need to use += operator (it is what works okay!). This modifies self's data. 
            # To ensure that the original affine expression's data is unmodified, temporary store data to save it.
            # This has to do with cython and overloaded operators. If this is issue is ever fixed, this could be revisited.  
            result_expr.thisptr = new Affine_Expr(self.thisptr[0]) 
            self.thisptr = temp
            return result_expr
        cdef Linear_Expr* e1
        if isinstance(other, Linear_Expression):
            other_lin = <Linear_Expression> other
            e1 = other_lin.thisptr
            self.thisptr[0] += e1[0]
            result_expr.thisptr = new Affine_Expr(self.thisptr[0])
            return result_expr
        # TODO: Make adding int explicit. 
        if isinstance(other, Affine_Expression):
            rhs = (<Affine_Expression> other).thisptr
        else:
            rhs_expr = Affine_Expression(other)
            rhs = (<Affine_Expression> rhs_expr).thisptr
        cdef Affine_Expr result
        result = lhs[0] + rhs[0]
        result_expr = Affine_Expression()
        result_expr.thisptr[0] = result
        return result_expr

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        r"""
        Subtract ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Affine_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Affine_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Affine_Expression`

        Examples:

        >>> from pplite import Affine_Expression, Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = x + y + y + y + 1; e
        x0+3*x1+1
        >>> e - e
        0
        >>> isinstance(e, Affine_Expression)
        True
        >>> 12 - e
        -x0-3*x1+11
        >>> e - 3 * Variable(3) - 3 * y
        x0-3*x3+1
        """
        cdef Affine_Expr* lhs
        if isinstance(self, Affine_Expression):
            lhs = (<Affine_Expression> self).thisptr
        else:
            lhs_expr = Affine_Expression(self)
            lhs = (<Affine_Expression> lhs_expr).thisptr
        result_expr = Affine_Expression()
        cdef Var* vv
        if isinstance(other, Variable):
            other_var = <Variable> other
            vv = other_var.thisptr
            self.thisptr[0] -= vv[0]
            result_expr.thisptr = new Affine_Expr(self.thisptr[0])
            return result_expr
        cdef Linear_Expr* e1
        if isinstance(other, Linear_Expression):
            other_lin = <Linear_Expression> other
            e1 = other_lin.thisptr
            self.thisptr[0] -= e1[0]
            result_expr.thisptr = new Affine_Expr(self.thisptr[0])
            return result_expr
        # TODO: Make adding int explicit.  
        if isinstance(other, Affine_Expression):
            rhs = (<Affine_Expression> other).thisptr
        else:
            rhs_expr = Affine_Expression(other)
            rhs = (<Affine_Expression> rhs_expr).thisptr
        cdef Affine_Expr result
        result = lhs[0] - rhs[0]
        result_expr = Affine_Expression()
        result_expr.thisptr[0] = result
        return result_expr

    def __rsub__(self, other):
        cdef Affine_Expr* lhs
        if isinstance(self, Affine_Expression):
            lhs = (<Affine_Expression> self).thisptr
        else:
            lhs_expr = Affine_Expression(self)
            lhs = (<Affine_Expression> lhs_expr).thisptr
        result_expr = Affine_Expression()
        cdef Var* vv
        if isinstance(other, Variable):
            other_var = <Variable> other
            vv = other_var.thisptr
            result_expr.thisptr = new Affine_Expr(vv[0] - self.thisptr[0]) 
            return result_expr
        # TODO: Make adding int explicit.
        if isinstance(other, Affine_Expression):
            rhs = (<Affine_Expression> other).thisptr
        else:
            rhs_expr = Affine_Expression(other)
            rhs = (<Affine_Expression> rhs_expr).thisptr
        cdef Affine_Expr result
        result = rhs[0] - lhs[0]
        result_expr = Affine_Expression()
        result_expr.thisptr[0] = result
        return result_expr

    def __mul__(self, other):
        cdef Affine_Expr* a
        if isinstance(self, Linear_Expression):
            a = (<Affine_Expression> self).thisptr
            c = other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            a = (<Affine_Expression> other).thisptr
            c = self
        cdef FLINT_Integer cc = Python_int_to_FLINT_Integer(c)
        cdef Affine_Expression result = Linear_Expression()
        result.thisptr[0] = a[0] * cc
        return result

    def __rmul__(self, other):
        cdef Affine_Expr* a
        if isinstance(self, Linear_Expression):
            a = (<Affine_Expression> self).thisptr
            c = other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            a = (<Affine_Expression> other).thisptr
            c = self
        cdef FLINT_Integer cc = Python_int_to_FLINT_Integer(c)
        cdef Affine_Expression result = Linear_Expression()
        result.thisptr[0] = a[0] * cc
        return result

    def __neg__(self):
        neg_assign(self.thisptr[0])
        return self

    def __pos__(self):
        return self

    def __richcmp__(self, other, op):
        """
        Construct :class:`Constraint`s

        Examples:

        >>> from pplite import Variable
        """
        return _make_Constraint_from_richcmp(self, other, op)