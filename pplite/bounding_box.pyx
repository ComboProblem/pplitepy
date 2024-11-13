
# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
# from libcpp.vector cimport vector as cppvector

from .integer_conversions cimport FLINT_Integer_to_Python, Python_int_to_FLINT_Integer, FLINT_Rational_to_Python, Python_float_to_FLINT_Rational
from .linear_algebra cimport Variable, Linear_Expression
from .constraint cimport Constraint
from .generators cimport PPliteGenerator
from .intervals cimport Interval

import_gmpy2()

cdef class Bounding_Box_t(object):
    """
    Wrapper for PPlite's BBox<true> class.

    Input: 

    Output: 
    
    TESTS:
    >>> from pplite.bounding_box import Bounding_Box_t
    >>> box_ex = Bounding_Box_t(1, "empty")
    >>> box_ex
    empty
    >>> box_ex_2 = Bounding_Box_t(box_ex); box_ex_2
    empty
    >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
    >>> x = Variable(0)
    >>> g = PPliteGenerator('ray', Linear_Expression(x), 1)
    >>> g_2 = PPliteGenerator('point', Linear_Expression(x), 0)
    >>> box_ex_2.add_generators([g,g_2]); box_ex_2 # Incorrect output, unsure why. might have something to do with interval class, conjecture either I"m not reading what should 
    # be defined here correctly. 
    x0 in [mpq(1,1), mpq(1,0)]
    # box isn't directly tested. it is tested indirectly via poly class. I should drop these tests. 

    """
    def __init__(self, *args): # special constructors as a class method. 
        if len(args) == 1:
            if isinstance(args[0], Bounding_Box_t):
                y = <Bounding_Box_t> args[0]
                self.thisptr = new Box_t(y.thisptr[0])
                return
            # if isinstance(other_box, Bounding_Box_f):
            #     y = <Bounding_Box_f> other_box
            #     self.thisptr = new Box_t(y.thisptr[0])
            #     return
        cdef dim_type sd
        cdef Spec_Elem se
        if len(args) == 2:
            sd = args[0]
            if args[1] == "empty":
                se = Spec_Elem.EMPTY
            elif args[1] == "universe":
                se = Spec_Elem.UNIVERSE
            else:
                raise ValueError("universe or empty box should be specified in second argument")
            self.thisptr = new Box_t(sd, se)
            return 

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        del self.thisptr

    def __repr__(self):
        s = ""
        if self.is_empty():
            s += "empty"
            return s
        else:
            for i in range(self.space_dimemsion()):
                s += str(Variable(i)) + " in " + str(self.get_bounds(Variable(i)))
                if i < self.space_dimemsion():
                    s += "\n"
            return s

    def __hash__(self):
        return self.thisptr.hash()

    def compute_volume_info(self):
        """
        Computes the volume of the Bounding Box.
        """
        cdef pair[dim_type, FLINT_Rational] Volume_Info
        Volume_Info = self.thisptr.compute_volume_info()
        cdef dim_type dim
        cdef FLINT_Rational vol
        dim = Volume_Info.first
        vol = Volume_Info.second
        py_vol = FLINT_Rational_to_Python(vol)
        return dim, py_vol

    def maybe_update_volume_info(self):
        self.thisptr.maybe_update_volume_info()

    def num_rays(self):
        """
        Returns the number of rays in the bounding box. 

        """
        return self.thisptr.num_rays()

    def pseudo_volume(self):
        """
        Returns the pseduo volume of the bounding box.
        """
        cdef FLINT_Rational p_vol
        p_vol = self.thisptr.pseudo_volume()
        return FLINT_Rational_to_Python(p_vol)

    def check_inv(self):
        return self.thisptr.check_inv()

    def space_dimemsion(self):
        return self.thisptr.space_dim()

    def affine_dimemsion(self):
        return self.thisptr.affine_dim()

    def is_empty(self):
        return self.thisptr.is_empty()

    def is_universe(self):
        return self.thisptr.is_universe()

    def is_bounded(self):
        return self.thisptr.is_bounded()

    def topology(self):
        pass

    def is_topologically_closed(self):
        return self.thisptr.is_topologically_closed()

    def constrains(self, variable):
        if isinstance(variable, Variable):
            vv = <Variable> variable
            return self.thisptr.constrains(vv.thisptr[0])

    def contains(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.contains(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.contains(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def equals(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.equals(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.equals(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def is_disjoint_from(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.is_disjoint_from(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.is_disjoint_from(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def less(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.less(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.less(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def number_of_minumum_constraints(self):
        return self.thisptr.num_min_cons()

    def number_of_minumum_generators(self):
        return self.thisptr.num_min_gens()

    def swap(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.swap(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.swap(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def set_empty(self):
        self.thisptr.set_empty()

    def set_origin(self):
        self.thisptr.set_origin()

    def affine_image(self, variable, expression, inhomogenous_term, denominator):
        cdef Var* var_ptr
        cdef Linear_Expr expr
        cdef FLINT_Integer inhomo
        cdef FLINT_Integer den
        # I want to use a try/except block but I havne't see code like that in sage and am unsure for this type ot type checking 
        # if the assumption that things are input correctly outweights teh patter of explicit checking. 
        if isinstance(variable, Variable):
            var_ptr = (<Variable> variable).thisptr
            if isinstance(expression, Linear_Expression):
                expr = (<Linear_Expression> expression).thisptr[0]
                inhom = Python_int_to_FLINT_Integer(inhomogenous_term)
                den = Python_int_to_FLINT_Integer(denominator)
                self.thisptr.affine_image(var_ptr[0], expr, inhomo, den)

    def concatenate_assign(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.concatenate_assign(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.concatenate_assign(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def glb_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.glb_assign(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def lub_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.lub_assign(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def unconstrain(self, variables):
        """
        INPUT: python iterable of :class:`Variable`
        """
        cdef Var* var
        for variable in variables:
            if isinstance(variable, Variable):
                var = (<Variable> variable).thisptr
                self.thisptr.unconstrain(var[0])

    def widening_assign(self, other_box):
        # if isinstance(other_box, Bounding_Box_f):
        #     w_o_vol_box = <Bounding_Box_f> other_box
        #     return self.thisptr.widening_assign(w_o_vol_box.thisptr[0])
        if isinstance(other_box, Bounding_Box_t):
            w_vol_box = <Bounding_Box_t> other_box
            return self.thisptr.widening_assign(w_vol_box.thisptr[0])

        raise ValueError("other_box needs to be a Bounding Box")

    def add_space_dims(self, dim, project):
        """
        Input: dim - int, project - bool
        """
        cdef dim_type d
        d = dim
        self.thisptr.add_space_dims(d, project)
        
    def remove_space_dim(self, variable):
        cdef Var* var
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
            self.thisptr.remove_space_dim(var[0])

    def remove_space_dims(self, variables):
        for variable in variables:
            self.remove_space_dim(variable)

    def remove_higher_space_dims(self, new_dimemsnion):
        cdef dim_type new_dim
        new_dim = new_dimemsnion
        self.thisptr.remove_higher_space_dims(new_dim)

    def expand_space_dim(self, variable, dimemsion_type):
        cdef Var* var
        cdef dim_type d
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
            d = dimemsion_type
            self.thisptr.expand_space_dim(var[0], d)

    def add_generator(self, generator):
        """
        """
        cdef Gen g
        if isinstance(generator, PPliteGenerator):
            g = (<PPliteGenerator> generator).thisptr[0]
            self.thisptr.add_gen(g)

    def add_generators(self, generators):
        for generator in generators:
            self.add_generator(generator)

    def refine_bounds(self, dimension, interval):
        cdef Itv y_itv
        if isinstance(interval, Interval):
            y_itv = (<Interval> interval).interval
            self.thisptr.refine_bounds(dimension, y_itv)

    def refine_as_intergral(self, *args):
        cdef dim_type dim
        cdef dim_type first
        cdef dim_type last
        if len(args) == 1:
            dim =  args[0]
            self.thisptr.refine_as_integral(dim)
        if len(args) == 2:
            first =  args[0]
            last = args[1]
            self.thisptr.refine_as_integral(first, last)            
        if len(args)>2:
            raise ValueError("expects 1 or 2 dimensions as input")

    def get_bounds(self, v):
        """
        """
        cdef Itv i
        if isinstance(v, Variable):
            vv = (<Variable> v).thisptr
            i = self.thisptr.get_bounds(vv[0])
            ii = Interval()
            ii.interval = i
            return ii

    def is_included_in(self, constraint):
        if isinstance(constraint, Constraint):
            cc = <Constraint> constraint
            return self.thisptr.is_included_in(cc.thisptr[0])

    def unconstrain_lb(self, dim):
        self.thisptr.unconstrain_lb(dim)

    def unconstrain_ub(self, dim):
        self.thisptr.unconstrain_ub(dim)

    def inf_lb(self, dim):
        return self.thisptr.inf_lb(dim)

    def inf_ub(self, dim):
        return self.thisptr.inf_ub(dim)

    def lb(self, dim):
        cdef FLINT_Rational q
        q = self.thisptr.lb(dim)
        return FLINT_Rational_to_Python(q)

    def ub(self, dim):
        cdef FLINT_Rational q
        q = self.thisptr.ub(dim)
        return FLINT_Rational_to_Python(q)

cdef class Bounding_Box_f(object):
    """
    Wrapper for PPlite's BBox<false> class.

    Input: 

    Output: 
    
    TESTS:
    >>> from pplite.bounding_box import Bounding_Box_f
    >>> box_ex = Bounding_Box_f(1, "empty")
    >>> box_ex
    empty
    >>> box_ex_2 = Bounding_Box_f(box_ex); box_ex_2
    empty
    >>> from pplite import Variable, Linear_Expression, Affine_Expression, PPliteGenerator
    >>> x = Variable(0)
    >>> g = PPliteGenerator('ray', Linear_Expression(x), 1)
    >>> g_2 = PPliteGenerator('point', Linear_Expression(x), 0)
    >>> box_ex_2.add_generators([g,g_2]); box_ex_2 # Incorrect output, unsure why. might have something to do with interval class, conjecture either I"m not reading what should 
    # be defined here correctly. 
    x0 in [mpq(1,1), mpq(1,0)]
    # box isn't directly tested. it is tested indirectly via poly class. I should drop these tests. 

    """
    def __init__(self, *args): # special constructors as a class method. 
        if len(args) == 1:
            # if isinstance(args[0], Bounding_Box_t):
            #     y = <Bounding_Box_t> args[0]
            #     self.thisptr = new Box_t(y.thisptr[0])
            #     return
            if isinstance(args[0], Bounding_Box_f):
                y = <Bounding_Box_f> args[0]
                self.thisptr = new Box_f(y.thisptr[0])
                return
        cdef dim_type sd
        cdef Spec_Elem se
        if len(args) == 2:
            sd = args[0]
            if args[1] == "empty":
                se = Spec_Elem.EMPTY
            elif args[1] == "universe":
                se = Spec_Elem.UNIVERSE
            else:
                raise ValueError("universe or empty box should be specified in second argument")
            self.thisptr = new Box_f(sd, se)
            return 

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        del self.thisptr

    def __repr__(self):
        s = ""
        if self.is_empty():
            s += "empty"
            return s
        else:
            for i in range(self.space_dimemsion()):
                s += str(Variable(i)) + " in " + str(self.get_bounds(Variable(i)))
                if i < self.space_dimemsion():
                    s += "\n"
            return s

    def __hash__(self):
        return self.thisptr.hash()

    def compute_volume_info(self):
        """
        Computes the volume of the Bounding Box.
        """
        cdef pair[dim_type, FLINT_Rational] Volume_Info
        Volume_Info = self.thisptr.compute_volume_info()
        cdef dim_type dim
        cdef FLINT_Rational vol
        dim = Volume_Info.first
        vol = Volume_Info.second
        py_vol = FLINT_Rational_to_Python(vol)
        return dim, py_vol

    def maybe_update_volume_info(self):
        self.thisptr.maybe_update_volume_info()

    def num_rays(self):
        """
        Returns the number of rays in the bounding box. 

        """
        return self.thisptr.num_rays()

    def pseudo_volume(self):
        """
        Returns the pseduo volume of the bounding box.
        """
        cdef FLINT_Rational p_vol
        p_vol = self.thisptr.pseudo_volume()
        return FLINT_Rational_to_Python(p_vol)

    def check_inv(self):
        return self.thisptr.check_inv()

    def space_dimemsion(self):
        return self.thisptr.space_dim()

    def affine_dimemsion(self):
        return self.thisptr.affine_dim()

    def is_empty(self):
        return self.thisptr.is_empty()

    def is_universe(self):
        return self.thisptr.is_universe()

    def is_bounded(self):
        return self.thisptr.is_bounded()

    def topology(self):
        pass

    def is_topologically_closed(self):
        return self.thisptr.is_topologically_closed()

    def constrains(self, variable):
        if isinstance(variable, Variable):
            vv = <Variable> variable
            return self.thisptr.constrains(vv.thisptr[0])

    def contains(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.contains(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.contains(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def equals(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.equals(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.equals(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def is_disjoint_from(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.is_disjoint_from(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.is_disjoint_from(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def less(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.less(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.less(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def number_of_minumum_constraints(self):
        return self.thisptr.num_min_cons()

    def number_of_minumum_generators(self):
        return self.thisptr.num_min_gens()

    def swap(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.swap(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.swap(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def set_empty(self):
        self.thisptr.set_empty()

    def set_origin(self):
        self.thisptr.set_origin()

    def affine_image(self, variable, expression, inhomogenous_term, denominator):
        cdef Var* var_ptr
        cdef Linear_Expr expr
        cdef FLINT_Integer inhomo
        cdef FLINT_Integer den
        # I want to use a try/except block but I havne't see code like that in sage and am unsure for this type ot type checking 
        # if the assumption that things are input correctly outweights teh patter of explicit checking. 
        if isinstance(variable, Variable):
            var_ptr = (<Variable> variable).thisptr
            if isinstance(expression, Linear_Expression):
                expr = (<Linear_Expression> expression).thisptr[0]
                inhom = Python_int_to_FLINT_Integer(inhomogenous_term)
                den = Python_int_to_FLINT_Integer(denominator)
                self.thisptr.affine_image(var_ptr[0], expr, inhomo, den)

    def concatenate_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.concatenate_assign(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.concatenate_assign(w_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def glb_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.glb_assign(w_o_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def lub_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.lub_assign(w_o_vol_box.thisptr[0])
        raise ValueError("other_box needs to be a Bounding Box")

    def unconstrain(self, variables):
        """
        INPUT: python iterable of :class:`Variable`
        """
        cdef Var* var
        for variable in variables:
            if isinstance(variable, Variable):
                var = (<Variable> variable).thisptr
                self.thisptr.unconstrain(var[0])

    def widening_assign(self, other_box):
        if isinstance(other_box, Bounding_Box_f):
            w_o_vol_box = <Bounding_Box_f> other_box
            return self.thisptr.widening_assign(w_o_vol_box.thisptr[0])
        # if isinstance(other_box, Bounding_Box_t):
        #     w_vol_box = <Bounding_Box_t> other_box
        #     return self.thisptr.widening_assign(w_vol_box.thisptr[0])

        raise ValueError("other_box needs to be a Bounding Box")

    def add_space_dims(self, dim, project):
        """
        Input: dim - int, project - bool
        """
        cdef dim_type d
        d = dim
        self.thisptr.add_space_dims(d, project)
        
    def remove_space_dim(self, variable):
        cdef Var* var
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
            self.thisptr.remove_space_dim(var[0])

    def remove_space_dims(self, variables):
        for variable in variables:
            self.remove_space_dim(variable)

    def remove_higher_space_dims(self, new_dimemsnion):
        cdef dim_type new_dim
        new_dim = new_dimemsnion
        self.thisptr.remove_higher_space_dims(new_dim)

    def expand_space_dim(self, variable, dimemsion_type):
        cdef Var* var
        cdef dim_type d
        if isinstance(variable, Variable):
            var = (<Variable> variable).thisptr
            d = dimemsion_type
            self.thisptr.expand_space_dim(var[0], d)

    def add_generator(self, generator):
        """
        """
        cdef Gen g
        if isinstance(generator, PPliteGenerator):
            g = (<PPliteGenerator> generator).thisptr[0]
            self.thisptr.add_gen(g)

    def add_generators(self, generators):
        for generator in generators:
            self.add_generator(generator)

    def refine_bounds(self, dimension, interval):
        cdef Itv y_itv
        if isinstance(interval, Interval):
            y_itv = (<Interval> interval).interval
            self.thisptr.refine_bounds(dimension, y_itv)

    def reinfe_as_integral(self, *args):
        cdef dim_type dim
        cdef dim_type first
        cdef dim_type last
        if len(args) == 1:
            dim =  args[0]
            self.thisptr.refine_as_integral(dim)
        if len(args) == 2:
            first =  args[0]
            last = args[1]
            self.thisptr.refine_as_integral(first, last)            
        if len(args)>2:
            raise ValueError("expects 1 or 2 dimensions as input")

    def get_bounds(self, v):
        """
        """
        cdef Itv i
        if isinstance(v, Variable):
            vv = (<Variable> v).thisptr
            i = self.thisptr.get_bounds(vv[0])
            ii = Interval()
            ii.interval = i
            return ii

    def is_included_in(self, constraint):
        if isinstance(constraint, Constraint):
            cc = <Constraint> constraint
            return self.thisptr.is_included_in(cc.thisptr[0])

    def unconstrain_lb(self, dim):
        self.thisptr.unconstrain_lb(dim)

    def unconstrain_ub(self, dim):
        self.thisptr.unconstrain_ub(dim)

    def inf_lb(self, dim):
        return self.thisptr.inf_lb(dim)

    def inf_ub(self, dim):
        return self.thisptr.inf_ub(dim)

    def lb(self, dim):
        cdef FLINT_Rational q
        q = self.thisptr.lb(dim)
        return FLINT_Rational_to_Python(q)

    def ub(self, dim):
        cdef FLINT_Rational q
        q = self.thisptr.ub(dim)
        return FLINT_Rational_to_Python(q)


# cpdef Create_Bounding_Box(keep_volume_info, arg_0 = None, arg_1 = None):
#     """
#     Wrapper for PPlite's BBox class. 

#     Input: Book keep_volue_info, 

#     Output: 

#     TESTS:

#     """
#     if arg_0 is None and arg_1 is None:
#         raise ValueError("At ")
#     if keep_volume_info == True:
#         return Bounding_Box_t(arg_0, arg_1)
#     if keep_volume_info == False:
#         return Bounding_Box_f(arg_0, arg_1)
# # #  