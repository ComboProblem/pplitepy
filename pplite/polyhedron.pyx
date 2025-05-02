# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
from libcpp.vector cimport vector as cppvector
from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer, FLINT_Rational_to_Python, Python_float_to_FLINT_Rational
from .constraint cimport Constraint
from .generators cimport PPliteGenerator
from .linear_algebra cimport Variable, Linear_Expression, Affine_Expression
from .intervals cimport Interval

cdef class NNC_Polyhedron(object):
    r"""
    Wrapper for PPLite's ``Poly`` class.
    
    INPUT:

    - dim_type (int), spec_elem (string "universe" or "empty"), and topology ("closed" or "nnc") xor,

    - nnc_poly, :class:`NNC_Polyhedron` xor,

    - cons - a list of :class:`Constraint`.

    OUTPUT:

    A :class:`NNC_Polyhedron`

    EXAMPLES::

    Construct an empty polyhedron:

    >>> from pplite import *
    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "empty", topology = "nnc")
    >>> P
    false

    Construct the universe:

    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
    >>> P

    Define the open first quadrant R^2 by adding constraints:

    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> P.add_constraint(x>0)
    >>> P.add_constraint(y>0)
    >>> P
    x0>0, x1>0

    We can define a polyhedron with another polyhedron:

    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> P.add_constraint(x>0)
    >>> P.add_constraint(y>0)
    >>> P
    x0>0, x1>0   
    >>> P_2 = NNC_Polyhedron(nnc_poly = P)
    >>> P_2
    x0>0, x1>0

    We can add generators in the form of points, closure points, lines, and rays to a polyhedron:

    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> P.add_constraint(x>0)
    >>> P.add_constraint(y>0)
    >>> P_2 = NNC_Polyhedron(nnc_poly = P)
    >>> P_2.add_generator(Point())  # This is a short cut for the origin
    >>> P_2
    x0>=0, x1>=0

    Add several constraints or generators at a time:

    >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> P.add_constraints([x>0, y>0])
    >>> P
    x0>0, x1>0

    Directly define a polyhedron from a list of constraints:

    >>> P_2 = NNC_Polyhedron(cons = [x>0, y>0])
    >>> P == P_2
    True    
    """
    def __init__(self, **kwrds):
        """
        TESTS::
        >>> from pplite import NNC_Polyhedron
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "empty", topology = "nnc")
        >>> P
        false
        >>> P_2 = NNC_Polyhedron(nnc_poly = P)
        >>> P_2
        false
        >>> P_2.equals(P)
        True
        >>> from pplite import Variable, Linear_Expression, Affine_Expression, Constraint
        >>> A = Variable(0)
        >>> P_3 = NNC_Polyhedron(dim_type = 1, spec_elem = "universe", topology = "nnc")
        >>> P_3.add_constraint(A >= 0)
        >>> P_3
        x0>=0
        >>> P_4 = NNC_Polyhedron(nnc_poly = P_3)
        >>> P_4.add_constraint(A >= -10)
        >>> P_4.add_constraint(A >= -5)
        >>> P_4.add_constraint(A >= -2)
        >>> P_4.add_constraint(A >= -1)
        >>> P_4.add_constraint(A >= 0)
        >>> P_3.equals(P_4)
        True
        >>> B = Variable(1)
        >>> cons_list = [A >= 0, B == 5]
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P.add_constraints(cons_list)
        >>> P
        x1-5==0, x0>=0
        """
        cdef dim_type dd
        cdef Spec_Elem ss
        cdef Topol tt
        cdef Poly* yy
        if len(kwrds) == 3:
            d = kwrds.pop("dim_type")
            s = kwrds.pop("spec_elem")
            t = kwrds.pop("topology")
            if isinstance(d, int):
                dd = d # needs to be a python int? Gotta figure out how to typed integer conversions. Ask about this. 
                ss = string_to_Spec_Elem(s)
                tt = string_to_Topol(t)
                self.thisptr = new Poly(dd, ss, tt)
                return
            raise ValueError("dim_type is required to be an instance of `:class:int`. spec_elem is either a string `universe` or `empty`, topology is either `nnc` or `closed`.")
        if "nnc_poly" in kwrds.keys():
            nnc_poly = kwrds.pop("nnc_poly")
            if isinstance(nnc_poly, NNC_Polyhedron):
                yy = (<NNC_Polyhedron> nnc_poly).thisptr
                self.thisptr = new Poly(yy[0])
                return
            raise TypeError(":class:`NNC_Polyhedron` needs to be provided to use the nnc_poly key word constructor.")
        if "cons" in kwrds.keys():
            cons = kwrds.pop("cons")
            d_cons = max([c.space_dimension() for c in cons])
            if "dim_type" in kwrds.keys():
                d_in = kwrds.pop("dim_type")
                d = max(d_cons, d_in)
            else:
                d = d_cons
            dd = d
            ss =  string_to_Spec_Elem("universe")
            
            tt = string_to_Topol("nnc")
            self.thisptr = new Poly(dd, ss, tt)
            for c in cons:
                cc = (<Constraint> c).thisptr[0] 
                self.thisptr.add_con(cc)
            return
        if "gens" in kwrds.keys():
            gens = kwrds.pop("gens")
            d_gens = max([g.space_dimension() for g in gens])
            if kwrds["dim_type"]:
                d_in = kwrds.pop("dim_type")
                d = max(d_gens, d_in)
            else:
                d = d_gens
            dd = d
            ss =  string_to_Spec_Elem("universe")
            tt = string_to_Topol("nnc")
            self.thisptr = new Poly(dd, ss, tt)
            for g in gens:
                gg = (<PPliteGenerator> g).thisptr[0] 
                self.thisptr.add_gen(gg)
            return
        raise KeyError("No valid keyword input for constructor of class given. ")

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        del self.thisptr

    def __hash__(self):
        return self.thisptr[0].hash()

    def __repr__(self):
        s = ""
        if self.is_empty(): 
            s = "false" # from pplite::Poly_Impl::print()
            return s
        self.minimize()
        comma = False
        for c in self.constraints():
            if comma:
                s += ", "
            s +=  str(c)
            comma = True
        return s

    def __eq__(self, other):
        if isinstance(other, NNC_Polyhedron):
            return self.equals(other)
        raise TypeError("Comparison with NNC polys only!")

    def is_necessarily_closed(self):
        return self.thisptr[0].is_necessarily_closed()

    def check_inv(self):
        return self.thisptr[0].check_inv()

    def is_empty(self):
        return self.thisptr[0].is_empty()

    def is_universe(self):
        return self.thisptr[0].is_universe()

    def is_minimized(self):
        return self.thisptr[0].is_minimized()

    def is_topologically_closed(self):
        return self.thisptr[0].is_topologically_closed()

    def is_bounded(self):
        return self.thisptr[0].is_bounded()

    def is_bounded_expression(self, from_below, expression):
        cdef Linear_Expr expr
        if isinstance(expression, Linear_Expression):
            expr = (<Linear_Expression> expression).thisptr[0]
        else:
            raise TypeError("expression needs to be of :class:`Linear_Expression`.")
        cdef cppbool f_b
        f_b = from_below
        return self.thisptr[0].is_bounded_expr(f_b, expr)

    def constrains(self,variable):
        cdef Var* vv
        if isinstance(variable, Variable):
            vv = (<Variable> variable).thisptr
            return self.thisptr[0].constrains(vv[0])
        raise TypeError("variable needs to be of :class:`Variable`.")

    def equals(self, other_poly):
        cdef Poly* yy
        if isinstance(other_poly, NNC_Polyhedron):
            yy = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].equals(yy[0])
        raise TypeError("other_poly needs to be of :class:`NNC_Polyhedron`")

    def contains(self, other_poly):
        cdef Poly* yy
        if isinstance(other_poly, NNC_Polyhedron):
            yy = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].contains(yy[0])
        raise TypeError("other_poly needs to be of :class:`NNC_Polyhedron`")

    def strictly_contains(self, other_poly):
        cdef Poly* yy
        if isinstance(other_poly, NNC_Polyhedron):
            yy = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].strictly_contains(yy[0])
        raise TypeError("other_poly needs to be of :class:`NNC_Polyhedron`")

    def is_disjoint_from(self, other_poly):
        cdef Poly* yy
        if isinstance(other_poly, NNC_Polyhedron):
            yy = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].is_disjoint_from(yy[0])
        raise TypeError("other_poly needs to be of :class:`NNC_Polyhedron`")

    def get_bounding_box(self):
        pass

    def boxed_contains(self, other_poly):
        cdef Poly* yy
        if isinstance(other_poly, NNC_Polyhedron):
            yy = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].boxed_contains(yy[0])
        raise TypeError("other_poly needs to be of :class:`NNC_Polyhedron`")

    # queries
    def topology(self):
        pass

    def space_dim(self):
        return self.thisptr.space_dim()

    def affine_dim(self):
        return self.thisptr[0].affine_dim()

    def _relation_with_c(self, constraint):
        cdef Con cc
        cdef Poly_Con_Rel p_c_r
        if isinstance(constraint, Constraint):
            cc = (<Constraint> constraint).thisptr[0]
            p_c_r = self.thisptr[0].relation_with(cc)
            result = Polyhedron_Constraint_Rel()
            result.thisptr[0] = p_c_r
            return result
        else:
            raise TypeError()

    def _relation_with_g(self, generator):
        cdef Gen gg
        cdef Poly_Gen_Rel p_g_r
        if isinstance(generator, PPliteGenerator):
            gg = (<PPliteGenerator> generator).thisptr[0]
            p_g_r = self.thisptr[0].relation_with(gg)
            result = Polyhedron_Generator_Rel()
            result.thisptr[0] = p_g_r
            return result
        else:
            raise TypeError("A :class:`PPliteGenerator` or a :class:`Constraint` should be passed into this method.")

    def relation_with(self, gen_or_constraint):
        try:
            return self._relation_with_c(gen_or_constraint)
        except TypeError:
            return self._relation_with_g(gen_or_constraint)  # failure here will raise the right type error for the general method. 



    def min(self, affine_expr, value, included_pointer, gen_object):
        cdef Affine_Expr ae
        cdef FLINT_Rational val
        cdef cppbool* included_ptr
        cdef Gen* g_ptr
        if isinstance(affine_expr, Affine_Expression):
            ae = (<Affine_Expression> affine_expr).thisptr[0]
        val = Python_float_to_FLINT_Rational(value)
        if isinstance(gen_object, PPliteGenerator):
            g_ptr = (<PPliteGenerator> gen_object).thisptr
        if isinstance(included_pointer, bool):
            included_ptr[0] = included_pointer # Doth this work?
        return self.thisptr[0].min(ae, val, included_ptr, g_ptr)

    def max(self, affine_expr, value, included_pointer, gen_object):
        cdef Affine_Expr ae
        cdef FLINT_Rational val
        cdef cppbool* included_ptr
        cdef Gen* g_ptr
        if isinstance(affine_expr, Affine_Expression):
            ae = (<Affine_Expression> affine_expr).thisptr[0]
        val = Python_float_to_FLINT_Rational(value)
        if isinstance(gen_object, PPliteGenerator):
            g_ptr = (<PPliteGenerator> gen_object).thisptr
        if isinstance(included_pointer, bool):
            included_ptr[0] = included_pointer # Doth this work?
        return self.thisptr[0].max(ae, val, included_ptr, g_ptr)

    def _get_bounds_v(self, variable):
        cdef Var* v
        cdef Itv itv 
        v = (<Variable> variable).thisptr
        itv = self.thisptr[0].get_bounds(v[0])
        i = Interval()
        i.interval = itv
        return i

    def _get_bounds_ae(self, affine_expr):
        cdef Affine_Expr ae
        cdef Itv itv 
        ae = (<Affine_Expression> affine_expr).thisptr[0]
        itv = self.thisptr[0].get_bounds(ae)
        i = Interval()
        i.interval = itv
        return i
        
    def _get_boundes_itv(self, itv_expr):
        pass

    def get_bounds(self, variable_or_affine_expr):
        if isinstance(variable_or_affine_expr, Variable):
            return self._get_bounds_v(variable_or_affine_expr)
        if isinstance(variable_or_affine_expr, Affine_Expression):
            return self._get_bounds_ae(variable_or_affine_expr)
        raise TypeError("A :class:`Variable` or a :class:`Affine_Expression` should be passed into this method.")

    def get_unconstrainted(self):
        pass

    def constraints(self):
        # Access constraints indirectly via copy_cons()
        # TODO: Properly implement via sys and Cons_Proxy in Poly_Impl
        cdef Cons constraint_vector 
        constraint_vector = self.thisptr[0].copy_cons()
        result = []
        cdef unsigned int index = constraint_vector.size() # hacky way to iterate over vectors
        for i in range(index):
            c = Constraint()
            c.thisptr = new Con(constraint_vector[i])
            result.append(c)
        return result

    def generators(self):
        """
        Returns a list of :class:`PPliteGenerator`.
        """
        # Access constraints indirectly via copy_cons()
        # TODO: Properly implement via sys and Cons_Proxy in Poly_Impl
        cdef Gens generator_vector 
        generator_vector = self.thisptr[0].copy_gens()
        result = []
        cdef unsigned int index = generator_vector.size()
        for i in range(index):
            g = PPliteGenerator()
            g.thisptr = new Gen(generator_vector[i])
            result.append(g)
        return result

    def normalized_constraints(self):
        # TODO implement once Cons_Proxy is implemented. 
        pass

    def num_min_constrains(self):
        return self.thisptr[0].num_min_cons()

    def num_min_generators(self):
        return self.thisptr[0].num_min_gens()

    def collapse(self, n):
        cdef dim_type nn
        if isinstance(n, int):
            nn = n
        self.thisptr[0].collapse(nn)

    def num_disjuncts(self):
        return self.thisptr[0].num_disjuncts()

    def disjunct_constraints(self, n):
        # TODO implement once Cons_Proxy is implemented. 
        # cdef dim_type nn
        # if isinstance(n, int):
        #     nn = n
        # cdef Cons_Proxy constraint_proxy
        pass

    def geom_covers(self, other_poly):
        """
        Input: :class:`NNC_Polyhedron`
        Output: bool
        """
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            return self.thisptr[0].geom_covers(y[0])

    def m_swap(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].m_swap(y[0])

    def set_empty(self):
        """
        Sets the :class:`NNC_Polyhedron` to empty. 
        """
        self.thisptr[0].set_empty()

    def set_topology(self, topology):
        cdef Topol tt
        tt = string_to_Topol(topology)
        self.thisptr[0].set_topology(tt)

    def add_constraint(self, constraint):  
        r"""
        TESTS::
        >>> from pplite import NNC_Polyhedron, Variable, Linear_Expression, Affine_Expression, Constraint, Point, Ray, Line, Closure_point
        >>> A = Variable(0)
        >>> P = NNC_Polyhedron(dim_type = 1, spec_elem = "universe", topology = "nnc")
        >>> P.add_constraint(A >= 0)
        >>> P
        x0>=0
        >>> P_2 = NNC_Polyhedron(nnc_poly = P)
        >>> P_2.add_constraint(A >= -10)
        >>> P_2.add_constraint(A >= -5)
        >>> P_2.add_constraint(A >= -2)
        >>> P_2.add_constraint(A >= -1)
        >>> P_2.add_constraint(A >= 0)
        >>> P_2.is_necessarily_closed()
        False
        >>> P_2.space_dim()
        1
        >>> P_2.equals(P)
        True
        >>> B = Variable(1)
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "empty", topology = "nnc") 
        >>> P_2 = NNC_Polyhedron(nnc_poly = P)
        >>> P.add_constraint(A == B)
        >>> P_2 == P
        True
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc") # addcons1 - test04
        >>> P.add_constraint(A >= 0)
        >>> P.add_constraint(A <= 2)
        >>> P.add_constraint(A >= -1)
        >>> P.add_constraint(B >= 1)
        >>> P_2 =  NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P_2.add_constraint(A >= 0)
        >>> P_2.add_constraint(A <= 2)
        >>> P_2.add_constraint(B >= 1)
        >>> P_2 == P
        True
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc") # addcons1 - test05
        >>> P.add_constraint(B >= 0)
        >>> P.add_constraint(B <= 2)
        >>> P.add_constraint(A + B >= 1)
        >>> P.add_constraint(A - B >= -1)
        >>> P_2 = NNC_Polyhedron(dim_type = 2, spec_elem = "empty", topology = "nnc")
        >>> P_2.add_generator(Point(B))
        >>> P_2.add_generator(Ray(A))
        >>> P_2.add_generator(Point(A+2*B))
        >>> P_2.add_generator(Point(A))
        >>> P == P_2
        True
        """
        if isinstance(constraint, Constraint):
            c = (<Constraint> constraint).thisptr[0]
            self.thisptr[0].add_con(c)
            return
        raise TypeError("Constraint not added. Check if input is a :class:`Constraint`.")


    def add_constraints(self, iter_of_cons):
        """
        TESTS::
        >>> from pplite import NNC_Polyhedron, Variable, Linear_Expression, Affine_Expression, Constraint
        >>> A = Variable(0)
        >>> B = Variable(1)
        >>> cons_list = [A >= 0, B == 5]
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc") #TODO: Make nicer python constructors
        >>> P.add_constraints(cons_list)
        >>> P_2 = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P_2.add_constraint(A >= 0)
        >>> P_2.add_constraint(B == 5)
        >>> P_2 == P
        True
        """
        for con in iter_of_cons:
            self.add_constraint(con)

    def add_generator(self, generator):
        """
        TESTS::
        >>> from pplite import NNC_Polyhedron, Variable, Linear_Expression, Affine_Expression, Constraint, Point, Closure_point, Ray, Line
        >>> A = Variable(0)
        >>> B = Variable(1)
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P.add_constraint(A >= 0)
        >>> P.add_constraint(B >= 0)
        >>> P.add_constraint(A + B > 0)
        >>> P.add_generator(Point())
        >>> P.minimize()
        >>> P_2 = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P_2.add_constraint(A >= 0)
        >>> P_2.add_constraint(B >= 0)
        >>> P_2.minimize()
        >>> P == P_2 # Test 01 finished
        True
        >>> P = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P.add_constraint(A >= 0)
        >>> P.add_constraint(B >= 0)
        >>> P.add_constraint(A + B > 0)
        >>> P.add_generator(Ray(-A))
        >>> P_2 = NNC_Polyhedron(dim_type = 2, spec_elem = "universe", topology = "nnc")
        >>> P_2.add_constraint(B >= 0)
        >>> P == P_2 # Test 03 finished
        True     
        """
        if isinstance(generator, PPliteGenerator):
            g = (<PPliteGenerator> generator).thisptr[0]
            self.thisptr[0].add_gen(g)

    def add_generators(self, iter_of_gens):
        for gen in iter_of_gens:
            self.add_generator(gen)

    def topological_closure_assign(self):
        self.thisptr[0].topological_closure_assign()

    def unconstain(self, variable):
        if isinstance(variable, Variable):
            v = (<Variable> variable).thisptr
            self.thisptr[0].unconstrain(v[0])

    def unconstain_many(self, iter_of_var_or_index_set):
        pass

    def intersection_assign(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].intersection_assign(y[0])

    def join_assign(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].join_assign(y[0])

    def poly_hull_assign(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].poly_hull_assign(y[0])

    def con_hull_assign(self, other_poly, boxed):
        cdef cppbool bboxed
        if boxed:
            bboxed = True
        else:
            bboxed = False
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].con_hull_assign(y[0], bboxed)

    def poly_difference_assign(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].poly_difference_assign(y[0])

    def affine_image(self, variable, linear_exp, inhomogenous_term, denominator):
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
        if isinstance(linear_exp, Linear_Expression):
            expr = (<Linear_Expression> linear_exp).thisptr[0]
        cdef FLINT_Integer inhomo
        cdef FLINT_Integer den
        inhomo = Python_int_to_FLINT_Integer(inhomogenous_term)
        den = Python_int_to_FLINT_Integer(denominator)
        self.thisptr[0].affine_image(var[0], expr, inhomo, den)

    def affine_preimage(self, variable, linear_exp, inhomogenous_term, denominator):
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
        if isinstance(linear_exp, Linear_Expression):
            expr = (<Linear_Expression> linear_exp).thisptr[0]
        cdef FLINT_Integer inhomo
        cdef FLINT_Integer den
        inhomo = Python_int_to_FLINT_Integer(inhomogenous_term)
        den = Python_int_to_FLINT_Integer(denominator)
        self.thisptr[0].affine_preimage(var[0], expr, inhomo, den)        
    
    # TODO: Implement these   
    def parallel_affine_image(self, args):
        pass

    def widing_assign(self, args):
        pass

    def time_elapse_assign(self, other_poly):
        if isinstance(other_poly, NNC_Polyhedron):
            y = (<NNC_Polyhedron> other_poly).thisptr
            self.thisptr[0].time_elapse_assign(y[0])

    def minimize(self):
        self.thisptr[0].minimize()

#####################################
### Poly_Con_Rel and Poly_Gen_Rel ###
#####################################

# TODO: Add full functionality of these classes.

cdef class Polyhedron_Constraint_Rel(object):
    def __cinit__(self):
        self.thisptr = NULL
    def __dealloc__(self):
        del self.thisptr

cdef class Polyhedron_Generator_Rel(object):
    def __cinit__(self):
        self.thisptr = NULL
    def __dealloc__(self):
        del self.thisptr


# TODO Migrate helper functions to a helper function module.

#########################
###  Helper Functions ###
#########################

cdef Topol string_to_Topol(t):
    cdef Topol tt
    if t == "closed":
        tt = Topol.CLOSED
        return tt
    if  t == "nnc":
        tt = Topol.NNC
        return tt
    raise ValueError("Topology type conversion failed.")

cdef Spec_Elem string_to_Spec_Elem(s):
    cdef Spec_Elem ss
    if s == "empty":
        ss = Spec_Elem.EMPTY
        return ss
    if  s == "universe":
        ss = Spec_Elem.UNIVERSE
        return ss
    raise ValueError("Spec_Elem type conversion failed.")