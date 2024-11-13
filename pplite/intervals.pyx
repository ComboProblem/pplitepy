# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from libcpp.vector cimport vector as cppvector
from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer, FLINT_Rational_to_Python, Python_float_to_FLINT_Rational
from .constraint cimport Constraint
from .linear_algebra cimport Variable

cdef class Interval(object):
    """
    cython wrapper for pplite ``Itv`` struct. 
    This represent a topologically closed 1 dimensional interval.

    EXAMPLES:
    >>> from pplite.intervals import Interval
    >>> I = Interval()
    >>> I["lb"] = 3/2
    >>> I["ub"] = 5/2
    >>> I
    [mpq(3,2), mpq(5,2)]
    >>> J = Interval()
    >>> J
    (-inf, +inf)

    Intervals are mutable:

    >>> I["lb"] = 4/2
    >>> I
    [mpq(2,1), mpq(5,2)]
    >>> I["lb"] = 100
    >>> I # this is a bug
    [mpq(100,1), mpq(5,2)] 
    >>> I.is_empty() # something is wrong here. investigate
    False
    """
    def __hash__(self):
        cdef size_t h
        h = self.interval.hash()
        return h

    def __repr__(self):
        """
        TESTS:
        """
        s = ""
        if self.is_universe():
            s += "(-inf, +inf)"
            return s
        elif self.inf_lb() and self.has_ub():
            s += "(-inf, "+self["ub"].__repr__()+"]"
            return s
        elif self.inf_ub() and self.has_lb():
            s += "["+self["lb"].__repr__()+", +inf)"
            return s
        elif self.has_ub() and self.has_lb():
            s += "["+self["lb"].__repr__()+", "+self["ub"].__repr__()+"]"
            return s
        elif self.is_empty():
            s += "{ }"
            return s

    def __setitem__(self, member, value):
        """
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["lb"]
        mpq(3,2)
        >>> I["ub"] = 12
        >>> I["ub"]
        mpq(12,1)
        """
        if member == "lb":
            self.set_lower_bound(value)
        if member == "ub":
            self.set_upper_bound(value)
        if member == "singlton":
            self.set_singleton(value)


    def __getitem__(self, member):
        """
        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["lb"]
        mpq(3,2)
        >>> I["ub"] = 12
        >>> I["ub"]
        mpq(12,1)
        """
        if member == "lb":
            return self.get_lower_bound()
        if member == "ub":
            return self.get_upper_bound()

    def set_lower_bound(self, value):
        """
        Sets the lower bound of the interval. 

        INPUT: Value - Rational or value convertible to python rational. 

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2 #  indirect test
        >>> I["lb"]
        mpq(3,2)
        """
        cdef FLINT_Rational lower_bound
        lower_bound = Python_float_to_FLINT_Rational(value)
        self.interval.set_lb(lower_bound)

    def get_lower_bound(self):
        """
        Gets the lower bound of the interval. 

        OUTPUT: Rational

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2 #  indirect test
        >>> I["lb"]
        mpq(3,2)
        """
        cdef FLINT_Rational lower_bound
        lower_bound = self.interval.lb
        return FLINT_Rational_to_Python(lower_bound)

    def set_upper_bound(self, value):
        """
        Sets the upper bound of the interval. 

        INPUT: Value - Rational or value convertible to python rational. 

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["ub"] = 3/2 #  indirect test
        >>> I["ub"]
        mpq(3,2)
        """
        cdef FLINT_Rational upper_bound
        upper_bound = Python_float_to_FLINT_Rational(value)
        self.interval.set_ub(upper_bound)

    def get_upper_bound(self):
        """
        Gets the upper bound of the interval. 

        OUTPUT: Rational

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["ub"] = 3/2 #  indirect test
        >>> I["ub"]
        mpq(3,2)
        """
        cdef FLINT_Rational upper_bound
        upper_bound  = self.interval.ub
        return FLINT_Rational_to_Python(upper_bound)

    def check_inv(self):
        """
        Calls the check inv method in the interval struct.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> I.check_inv()
        True
        """
        return self.interval.check_inv()

    def is_empty(self):
        """
        Checks if the interval is empty.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> I.is_empty()
        False
        >>> J = Interval()
        >>> J.set_empty()
        >>> J.is_empty()
        True
        """
        return self.interval.is_empty()

    def is_universe(self):
        """
        Checks if the interval is universe (i.e. the real line).

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> I.is_universe()
        False
        >>> J = Interval()
        >>> J.set_universe()
        >>> J.is_universe()
        True
        """
        return self.interval.is_universe()

    def has_lb(self):
        """
        Checks if the interval has a lower bound.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I.has_lb()
        True
        >>> J = Interval()
        >>> J["ub"] = 12
        >>> J.has_lb()
        False
        """
        return self.interval.has_lb()

    def has_ub(self):
        """
        Checks if the interval has an upper bound.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I.has_ub()
        False
        >>> J = Interval()
        >>> J["ub"] = 12
        >>> J.has_ub()
        True
        """
        return self.interval.has_ub()

    def inf_lb(self):
        """
        Checks if the interval is unbounded from below.
        Output: Bool

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I.inf_lb()
        True
        """
        return self.interval.inf_lb()

    def inf_ub(self):
        """
        Checks if the interval is unbounded from above.
        Output: Bool

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I.inf_ub()
        True
        """
        return self.interval.inf_ub()

    def is_bounded(self):
        """
        Checks if the interval has an upper bound.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I.has_lb()
        True
        >>> J = Interval()
        >>> J["ub"] = 12
        >>> J.has_lb()
        False
        """
        return self.interval.is_bounded()

    def is_singleton(self):
        """
        Checks if the interval has an upper bound.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 3/2
        >>> I.is_singleton()
        True
        >>> J = Interval()
        >>> J["ub"] = 3/2
        >>> J["ub"] = 12
        >>> J.is_singleton()
        False
        """
        return self.interval.is_singleton()

    def is_zero(self):
        """
        Checks if the interval is zero.

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> I.is_zero()
        False
        >>> J = Interval()
        >>> J.set_zero()
        >>> J.is_zero()
        True
        """
        return self.interval.is_zero()

    def is_disjoint_from(self, other_interval):
        """
        Determines if interval I is disjoint form interval J. 

        INPUT: other_interval a :class:`Interval`

        OUTPUT: Bool
        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> J = Interval()
        >>> J.set_zero()
        >>> I.is_disjoint_from(J)
        True
        >>> J["ub"] = 10
        >>> I.is_disjoint_from(J)
        False
        """
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            return self.interval.is_disjoint_from(y)
        raise ValueError("other_interval is required to be an Interval")

    def intersects(self, other_interval):
        """
        Determines if interval I is intersects interval J. 

        INPUT: other_interval a :class:`Interval`

        OUTPUT: Bool

        TESTS: 
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 12
        >>> J = Interval()
        >>> J.set_zero()
        >>> I.intersects(J)
        False
        >>> J["ub"] = 10
        >>> I.intersects(J)
        True
        """
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            return self.interval.intersects(y)
        raise ValueError("other_interval is required to be an Interval")

    def length(self):
        """
        Computes the length of an interval.

        OUTPUT: :class:`mpq`

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I["lb"] = 3/2
        >>> I["ub"] = 5/2
        >>> I.length()
        mpq(1,1)
        """
        cdef FLINT_Rational l 
        l = self.interval.length()
        return FLINT_Rational_to_Python(l)

    def num_min_cons(self):
        """

        """
        cdef dim_type d
        d = self.interval.num_min_cons()
        return d

    def num_rays(self):
        """
        Returns the number of rays.
        """
        cdef dim_type d
        d = self.interval.num_rays()
        return d

    def set_empty(self):
        """
        Sets the interval to the empty set.

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval(); I
        (-inf, +inf)
        >>> I.set_empty()
        >>> I
        { }
        """
        self.interval.set_empty()

    def set_universe(self):
        """
        Sets the interval to the empty set.

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval();
        >>> I["ub"] = 12; I["lb"] = 2; I
        [mpq(2,1), mpq(12,1)]
        >>> I.set_universe()
        >>> I
        (-inf, +inf)
        """
        self.interval.set_universe()

    def set_zero(self):
        """
        Sets the interval to the empty set.

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I
        (-inf, +inf)
        >>> I.set_zero()
        >>> I
        [mpq(0,1), mpq(0,1)]
        """
        self.interval.set_zero()

    def set_singleton(self, value):
        """
        Sets the interval to a singleton of the specified value.

        INPUT: rational

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval()
        >>> I
        (-inf, +inf)
        >>> I.set_singleton(3/4)
        >>> I
        [mpq(3,4), mpq(3,4)]
        """
        cdef FLINT_Rational v
        v = Python_float_to_FLINT_Rational(value)
        self.interval.set_singleton(v)

    def unset_lower_bound(self):
        """
        Unset the interval's lower bound value.

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval();
        >>> I["ub"] = 12; I["lb"] = 2; I
        [mpq(2,1), mpq(12,1)]
        >>> I.unset_lower_bound()
        >>> I
        (-inf, mpq(12,1)]
        """
        self.interval.unset_lb()

    def unset_upper_bound(self):
        """
        Unset the interval's upperbound bound value.

        TESTS:
        >>> from pplite.intervals import Interval
        >>> I = Interval();
        >>> I["ub"] = 12; I["lb"] = 2; I
        [mpq(2,1), mpq(12,1)]
        >>> I.unset_upper_bound()
        >>> I
        [mpq(2,1), +inf)
        """
        self.interval.unset_ub()

# Double check what these functions should do. 

    def glb_assign(self, other_interval):
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            return self.interval.glb_assign(y)
        raise ValueError("other_interval must be an Interval.")

    def lub_assign(self, other_interval):
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            self.interval.lub_assign(y)
        raise ValueError("other_interval must be an Interval.")            

    def widen_assign(self, other_interval):
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            self.interval.widen_assign(y)
        raise ValueError("other_interval must be an Interval.")

    def refine_as_integral(self):
        return self.interval.refine_as_integral()

    def complement_assign(self):
        self.interval.complement_assign()

    def add_assign(self, other_interval):
        cdef Itv y
        if isinstance(other_interval, Interval):
            y = (<Interval> other_interval).interval
            self.interval.add_assign(y)
        raise ValueError("other_interval must be an Interval.")

    def mul_assign(self, value):
        cdef FLINT_Rational r
        r = Python_float_to_FLINT_Rational(value)
        self.interval.mul_assign(r)

###########################
###########################
### Functions in itv.hh ###
###########################
###########################

def interval_from_con_inhomo(constraint):
    cdef Con c
    cdef Itv result
    if isinstance(constraint, Constraint):
        c = (<Constraint> constraint).thisptr[0]
        result = itv_from_con_inhomo(c)
        new_interval = Interval()
        new_interval.interval = result
        return new_interval
    raise ValueError("constraint must be a Constraint.")

def interval_from_itv_con(constraint):
    cdef Con c
    cdef Itv result
    if isinstance(constraint, Constraint):
        c = (<Constraint> constraint).thisptr[0]
        result = itv_from_itv_con(c)
        new_interval = Interval()
        new_interval.interval = result
        return new_interval
    raise ValueError("constraint must be a Constraint.")

def split_interval(this_interval, constraint, intergral):
    """
    INPUT: this_interval :class:`Interval`, constraint :class:`Constraint`, integral :bool:
    OUTPUT: :class:`Interval`
    """
    cdef Itv itv
    cdef Con c
    cdef cppbool b
    cdef Itv result
    if isinstance(this_interval, Interval):
        itv = (<Interval> this_interval).interval
    if isinstance(constraint, Constraint):
        c = (<Constraint> constraint).thisptr[0]
    if intergral is True:
        b = True
    else:
        b = False
    result = split_itv(itv, c, b)
    new_interval = Interval()
    new_interval.interval = result
    return new_interval

def get_lower_bound_constraint(variable, this_interval):
    cdef Var* var
    cdef Itv itv
    cdef Con result
    if isinstance(variable, Variable):
        var = (<Variable> variable).thisptr
    else:
        var = (<Variable> Variable(variable)).thisptr
    if isinstance(this_interval, Interval):
        itv = (<Interval> this_interval).interval
    else:
        raise ValueError("this_interval is required to be an Interval")
    result = get_lb_con(var[0], itv)
    new_constraint = Constraint()
    new_constraint.thisptr[0] = result
    return new_constraint

def get_upper_bound_constraint(variable, this_interval):
    cdef Var* var
    cdef Itv itv
    cdef Con result
    if isinstance(variable, Variable):
        var = (<Variable> variable).thisptr
    else:
        var = (<Variable> Variable(variable)).thisptr
    if isinstance(this_interval, Interval):
        itv = (<Interval> this_interval).interval
    else:
        raise ValueError("this_interval is required to be an Interval")
    result = get_ub_con(var[0], itv)
    new_constraint = Constraint()
    new_constraint.thisptr[0] = result
    return new_constraint

def get_equality_constraint(variable, this_interval):
    cdef Var* var
    cdef Itv itv
    cdef Con result
    if isinstance(variable, Variable):
        var = (<Variable> variable).thisptr
    else:
        var = (<Variable> Variable(variable)).thisptr
    if isinstance(this_interval, Interval):
        itv = (<Interval> this_interval).interval
    else:
        raise ValueError("this_interval is required to be an Interval")
    result = get_eq_con(var[0], itv)
    new_constraint = Constraint()
    new_constraint.thisptr[0] = result
    return new_constraint