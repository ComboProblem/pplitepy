# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

from gmpy2 cimport GMPy_MPZ_From_mpz, import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from cython.operator cimport dereference as deref
from .linear_algebra import Variable, Linear_Expression, Affine_Expression
from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer

cdef class PPliteGenerator(object):
    """
    Wrapper for PPlite's Gen class.

    Input: :class:`PPliteGenerator`, or Generator Type, :class:`Linear_Experssion` and `int`.

    Output: :class:`PPlite_Generator`
    """
    def __init__(self, *args):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> PPliteGenerator()
        p(0)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('line', e, 7); g
        l(2*x0+3*x1)
        >>> g_2 = PPliteGenerator(g); g_2
        l(2*x0+3*x1)
        >>> e_2 = 2*x-3*y
        >>> g_3 =  PPliteGenerator('point', e_2, 9); g_3
        p((2*x0-3*x1)/9)
        >>> g_4 = PPliteGenerator('ray', Linear_Expression(x), 0); g_4
        r(x0)
        """
        if len(args) == 0:
            self.thisptr = new Gen()
            return
        cdef Gen gg
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, PPliteGenerator): 
                gg = (<PPliteGenerator> arg).thisptr[0]
                self.thisptr = new Gen(gg)
                return
            raise ValueError("creating an instance of PPliteGenerator from one input requires the input to be a PPliteGenerator.")
        cdef GenType tt
        cdef FLINT_Integer dd
        cdef Linear_Expr ee 
        if len(args) == 3:
            t = args[0]
            e = args[1]
            d = args[2]
            tt = string_to_GenType(t)
            if isinstance(e, Linear_Expression):
                ee = (<Linear_Expression> e).thisptr[0]
                dd = Python_int_to_FLINT_Integer(d)
                self.thisptr = new Gen(tt, ee, dd)
                return
            raise ValueError("Three inputs requires a Generator Type, linear expression, and integer")
        raise ValueError("Four or more arguments cannot be used to define a generator.")

    def __dealloc__(self):
        """
        The cython deconstructor.
        """
        del self.thisptr

    def __repr__(self):
        """
        Representation of a generator in the space of generators. 

        For example p(x0+x1) corresponds to the point (1,1) as an ordered pair. 

        OUTPUT::

        string
        """
        need_div = False
        s = ""
        if self.gen_type() == "line":
            s += "l("
        elif self.gen_type() == "ray":
            s += "r("
        elif self.gen_type() == "point":
            s += "p("
            if self.divisor() != 1:
                need_div = True
        elif self.gen_type() == "closure_point":
            s += "c("
            if self.divisor() != 1:
                need_div = True
        if need_div:
            s += "("
        s += str(self.linear_form())
        if need_div:
            s += ")/"
            s += str(self.divisor())
        s += ")"
        return s

    def __hash__(self):
        raise TypeError("PPliteGenerator unhashable")

    def __eq__(self, other_gen):
        cdef Gen y
        if isinstance(other_gen, PPliteGenerator):
            y =  (<PPliteGenerator> other_gen).thisptr[0]
            return self.thisptr[0] == y

    def gen_type(self):
        """
        Returns the type of generator. 

        Output: "line", "ray", "point", or "closure_point"

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('line', e, 7)
        >>> g.gen_type()
        'line'
        """
        cdef GenType t
        t = self.thisptr.type()
        if t == GenType.LINE:
            return "line"
        if t == GenType.RAY:
            return "ray"
        if t == GenType.POINT:
            return "point"
        if t == GenType.CLOSURE_POINT:
            return "closure_point"

    def set_type(self, t):
        """
        Sets the generator type of a generatorg

        Output: "line", "ray", "point", or "closure_point"

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('line', e, 7)
        >>> g.gen_type()
        'line'
        >>> g.set_type("ray")
        >>> g.gen_type()
        'ray'
        """
        cdef GenType tt
        tt = string_to_GenType(t)
        self.thisptr.set_type(tt)

    def is_line(self):
        """
        Returns true if the generator is a line. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('line', e, 7)
        >>> g.is_line()
        True
        >>> g = PPliteGenerator('ray', e, 7)
        >>> g.is_line()
        False
        """
        return self.thisptr.is_line()

    def is_ray(self):
        """
        Returns true if the generator is a ray. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('line', e, 7)
        >>> g.is_ray()
        False
        >>> g = PPliteGenerator('ray', e, 7)
        >>> g.is_ray()
        True
        """
        return self.thisptr.is_ray()

    def is_point(self):
        """
        Returns true if the generator is a point. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.is_point()
        True
        >>> g = PPliteGenerator('closure_point', e, 7)
        >>> g.is_point()
        False
        """
        return self.thisptr.is_point()

    def is_closure_point(self):
        """
        Returns true if the generator is a closure point. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.is_closure_point()
        False
        >>> g = PPliteGenerator('closure_point', e, 7)
        >>> g.is_closure_point()
        True
        """
        return self.thisptr.is_closure_point()

    def is_line_or_ray(self):
        """
        Returns true if the generator is a line or ray. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.is_line_or_ray()
        False
        >>> g = PPliteGenerator('line', e, 7)
        >>> g.is_line_or_ray()
        True
        >>> g = PPliteGenerator('ray', e, 7)
        >>> g.is_line_or_ray()
        True
        """
        return self.thisptr.is_line_or_ray()

    def is_point_or_closure_point(self):
        """
        Returns true if the generator is a point or closure point. 

        Output: bool

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.is_point_or_closure_point()
        True
        >>> g = PPliteGenerator('closure_point', e, 7)
        >>> g.is_point_or_closure_point()
        True
        >>> g = PPliteGenerator('ray', e, 7)
        >>> g.is_point_or_closure_point()
        False
        """
        return self.thisptr.is_point_or_closure_point()

    def linear_form(self):
        """
        Returns the Linear Form associated with the generator.

        Output: :class:`Linear_Expression`

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.linear_form()
        2*x0+3*x1
        """
        cdef Linear_Expr e
        e = self.thisptr.linear_expr()
        expr = Linear_Expression()
        expr.thisptr[0] = e
        return expr

    def coefficient(self, v):
        """
        Returns the coefficient of variable for a given generator. 

        Input: `int` or :class:`Variable`

        Output: `int`

        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.coefficient(x)
        mpz(2)
        """
        cdef Var* vv
        if isinstance(v, Variable):
            vv = (<Variable> v).thisptr
        else:
            vv = (<Variable> Variable(v)).thisptr
        cdef FLINT_Integer n
        n = self.thisptr.coeff(vv[0])
        return FLINT_Integer_to_Python(n)

    def divisor(self):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.divisor()
        mpz(7)
        """
        cdef FLINT_Integer n
        n = self.thisptr.divisor()
        return FLINT_Integer_to_Python(n)

    def space_dimision(self):
        """
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> g = PPliteGenerator('point', e, 7)
        >>> g.space_dimision()
        2
        """
        return self.thisptr.space_dim()

    def set_space_dimension(self, dim_type dim):
        self.thisptr.set_space_dim(dim)

    def is_equal_to(self, y):
        """
        Tests if self is equal to y.
        TESTS:
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x + 3*y
        >>> e_2 = x + y
        >>> g_1 = PPliteGenerator('point', e, 7)
        >>> g_2 = PPliteGenerator('ray', e, 7)
        >>> g_1.is_equal_to(g_2)
        False
        >>> g_3 = PPliteGenerator('point', e, 7)
        >>> g_1.is_equal_to(g_3)
        True
        >>> g_4 = PPliteGenerator('line', e_2, 4)
        >>> g_4.is_equal_to(g_3)
        False
        """
        cdef Gen gg
        if isinstance(y, PPliteGenerator):
            gg = (<PPliteGenerator> y).thisptr[0]
            return self.thisptr.is_equal_to(gg)
        raise ValueError("This is for comparing two generators.")

    def check_inv(self):
        return self.thisptr.check_inv()

    def m_swap(self, y):
        cdef Gen gg
        if isinstance(y, PPliteGenerator):
            gg = (<PPliteGenerator> y).thisptr[0]
            self.thisptr.m_swap(gg)
            return
        raise ValueError("This works only when using an input of a generator")

    def permute_space_dims(self, cycle, d):
        return NotImplemented

    def shift_space_dims(self, dim_type start, dim_type n):
        self.thisptr.shift_space_dims(start, n)

    def is_line_or_equality(self):
        return self.thisptr.is_line_or_equality()

    def sign_normalize(self):
        self.thisptr.sign_normalize()

    def strong_normalize(self):
        self.thisptr.strong_normalize()

###########################
### Auxiliary Functions ###
###########################

# Provides equivalent functionality to line, ray, point, closure_point functions in gen.hh

def Point(linear_expression=Linear_Expression(), integer=1):
    """
    Constructs a point.
    """
    return PPliteGenerator("point", Linear_Expression(linear_expression), integer)

def Closure_point(linear_expression=Linear_Expression(), integer=1):
    """
    Constructs a closure point.
    """
    return PPliteGenerator("closure_point", Linear_Expression(linear_expression), integer)

def Line(linear_expression):
    """
    Constructs a line.
    """
    return PPliteGenerator("line", Linear_Expression(linear_expression), 0)

def Ray(linear_expression):
    """
    Constructs a ray.
    """
    return PPliteGenerator("ray", Linear_Expression(linear_expression), 0)
    

########################
### Helper Functions ###
########################

# TODO migrate helper functions to a single module. 

cdef GenType string_to_GenType(t):
    """
    Converts a string to an enum pplite::GenType.

    Input: String - 'line', 'ray', 'point', 'closure point'
    
    Output: GenType of the same name.
    """
    cdef GenType tt
    if t == "line":
        tt = GenType.LINE
        return tt
    if t == "ray":
        tt = GenType.RAY
        return tt
    if t == "point":
        tt = GenType.POINT
        return tt
    if t == "closure_point":
        tt = GenType.CLOSURE_POINT
        return tt
    raise ValueError("Input must be one of the following string: `line`, `ray`, `point`, `closure point`.")