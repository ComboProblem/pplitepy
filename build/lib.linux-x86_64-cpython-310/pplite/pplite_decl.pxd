from libcpp cimport bool as cppbool
from libcpp.vector cimport vector as cppvector
from libcpp.utility cimport pair
from gmpy2 cimport mpz 

# gmp and flint integer/rational cdefs 

cdef extern from "gmp.h":
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_limb_signed_t
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr
    void mpz_init(mpz_t)
    void mpz_clear(mpz_t)
    cdef mpz_t* address_of_mpz "&"(mpz_t x)

cdef extern from "gmpxx.h":
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)

ctypedef mp_limb_t ulong
ctypedef mp_limb_signed_t slong

cdef extern from "flint/fmpz.h":
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]
    void fmpz_get_mpz(mpz_t x, const fmpz_t f) noexcept
    void fmpz_init(fmpz_t f)
    void fmpz_set_ui(fmpz_t f, ulong g)
    void fmpz_set_si(fmpz_t f, slong g)
    void fmpz_clear(fmpz_t f)
    int fmpz_print(const fmpz_t x)

cdef extern from "flint/fmpq.h":
    ctypedef struct fmpq:
        pass
    ctypedef fmpq fmpq_t[1]
    void fmpz_init(fmpq_t x)
    void fmpq_clear(fmpq_t x)
    void fmpq_set(fmpq_t dest, const fmpq_t src)
    void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)
    void fmpq_set_si(fmpq_t res, slong p, ulong q)

# Starting pplite definitons

cdef extern from "pplite/pplite.hh" namespace "pplite":

# "pplite/FLINT_Integer.hh":

    cdef cppclass FLINT_Integer
    cdef cppclass FLINT_Integer:
        FLINT_Integer()
        const fmpz_t& impl() 
        FLINT_Integer(const fmpz_t z)
        FLINT_Integer(signed int si)
        FLINT_Integer(const mpz_t z)

# "pplite/FLINT_Rational.hh":
    
    cdef cppclass FLINT_Rational
    cdef cppclass FLINT_Rational:
        FLINT_Rational()
        FLINT_Rational(FLINT_Rational x)
        FLINT_Rational(const FLINT_Integer& n, const FLINT_Integer& d)
        fmpq* impl() 
        void operator=(FLINT_Rational x)


# "pplite/GMP_Integer.hh"

#     cdef cppclass GMP_Integer
#     cdef cppclass GMP_Integer:
#         GMP_Integer()
#         const mpz_t& impl() 
#         # GMP_Integer(const fmpz_t z)
#         GMP_Integer(signed int si)
#         GMP_Integer(const mpz_t z)

# "pplite/globals.hh"
 
    ctypedef size_t dim_type
    cdef enum class Spec_Elem:
        EMPTY
        UNIVERSE

    cdef enum class Topol:
        CLOSED
        NNC

    # cdef enum class Widen_Spec:
    #     SAFE
    #     RISKY

    # cdef enum class Widen_Impl:
    #     H79
    #     BOXED_H79
    #     BHRZ03

    # Topol get_default_topology()


# "pplite/Var.hh"

    cdef cppclass Var
    cdef cppclass Var:
        Var(dim_type i)
        dim_type id()
        dim_type space_dim()
        void m_swap(Var& v)
    ctypedef struct Vars_Set
    cppbool less_than(Var v, Var w)
    void swap(Var& v, Var& w)
    
# "pplite/Linear_Expr.hh"

    ctypedef cppvector[FLINT_Integer] Impl # I might need to move this to the linear alg. pyx file
    cdef cppclass Linear_Expr
    cdef cppclass Linear_Expr:
        Linear_Expr()
        Linear_Expr(Var v)
        Linear_Expr(dim_type dim)
        Linear_Expr(Linear_Expr &e)
        Linear_Expr(const Linear_Expr &e, dim_type dim)
        # Linear_Expr operator=(Linear_Expr &e)
        Linear_Expr operator+(Linear_Expr &e)
        # Linear_Expr operator+(Var v, Linear_Expr e)
        Linear_Expr operator-(Linear_Expr &e)
        Linear_Expr operator*(FLINT_Integer &c)
        dim_type id()  #methods? #linear_expr.hh lines 39-112
        dim_type space_dim()
        FLINT_Integer get(dim_type dim)
        FLINT_Integer get(Var v)
        void set(dim_type dim, FLINT_Integer n)
        void set(Var v, FLINT_Integer n)
        void set_space_dim(dim_type dim)
        void swap_space_dims(dim_type i, dim_type j)
        void shift_space_dims(dim_type start, dim_type n)
        void shift_space_dims(Var start, dim_type n)
        # void remove_space_dims(Iter first, Iter last)
        Impl impl()
        cppbool is_equal_to(Linear_Expr& y) 
        cppbool is_zero()
    # Linear_Expr operator+(Linear_Expr e1, Linear_Expr e2)

    Linear_Expr operator+(Linear_Expr e, Var v) 
    Linear_Expr operator+(Var v, Linear_Expr v)
    Linear_Expr operator+(Var v, Var w)
    # Linear_Expr operator-(Linear_Expr e1, Linear_Expr e2)
    Linear_Expr operator-(Linear_Expr &e, Var &v)
    Linear_Expr operator-(Var &v, Linear_Expr &v)
    Linear_Expr operator-(Var v, Var w)
    void neg_assign(Linear_Expr& e)

# "pplite/Affine_Expr.hh"

    cdef cppclass Affine_Expr
    cdef cppclass Affine_Expr:
        Linear_Expr expr
        FLINT_Integer inhomo
        Affine_Expr()
        Affine_Expr(FLINT_Integer i)
        Affine_Expr(Linear_Expr e, FLINT_Integer i)
        Affine_Expr(Affine_Expr &e)
        dim_type space_dim()
        void set_space_dim(dim_type dim)
        cppbool is_zero()
        void m_swap(Affine_Expr& y)
        void normalize()
        void sign_normalize()
        Affine_Expr operator+(Affine_Expr &a)
    #    Affine_Expr operator+(Linear_Expr &e, FLINT_Integer &n)
        Affine_Expr operator-(Affine_Expr &a)
        Affine_Expr operator*(FLINT_Integer &c)
    # Affine_Expr operator+(Linear_Expr &e, Var &v)
    # Affine_Expr operator+(Linear_Expr &e, FLINT_Integer &n)
    # Affine_Expr operator+(Affine_Expr &e, FLINT_Integer &n)
    # Affine_Expr operator+(Affine_Expr &a1, Linear_Expr &e2)
    # Affine_Expr operator+(Affine_Expr &a, Linear_Expr &e2)
    Affine_Expr operator+(Affine_Expr a, Var v)
    Affine_Expr operator-(Var v, Affine_Expr a)
    Affine_Expr operator-(Linear_Expr e1, Affine_Expr a1)
    void neg_assign(Affine_Expr& a)
    # Affine_Expr& operator+=(Affine_Expr& a1, Var v) #+= operator not yet supported.
     
    # "pplite/Con.hh"

    cdef cppclass Con
    cdef cppclass Con:
        enum ConType "Con::Type":
            EQUALITY
            NONSTRICT_INEQUALITY
            STRICT_INEQUALITY
        struct Impl:
            Linear_Expr expr
            FLINT_Integer inhomo
            ConType type
            Impl()
            Impl(Linear_Expr e, FLINT_Integer i, ConType t)
        Con()
        Con(const Con &c)
        Con(Linear_Expr expr, FLINT_Integer inhomo, ConType type)
        Con(Affine_Expr ae, ConType type)
        dim_type space_dim()    
        void set_space_dim(dim_type dim)
        # void permute_space_dims_cycle(const Dims& cycle, dim_type d)
        void shift_space_dims(dim_type start, dim_type n) 
        Impl& impl()
        ConType type()
        cppbool is_equality()
        cppbool is_inequality()
        cppbool is_nonstrict_inequality()
        cppbool is_strict_inequality()
        Linear_Expr linear_expr()
        FLINT_Integer coeff(Var v)
        FLINT_Integer inhomo_term()
        Con zero_dim_false()
        Con zero_dim_positivity()
        cppbool is_tautological()
        cppbool is_inconsistent()
        cppbool is_equal_to(const Con& y)
        cppbool check_inv()
        void m_swap(Con& y)
        void set_type(ConType t)
        cppbool is_line_or_equality()
        void set_is_line_or_equality()
        void linear_combine(const Con& y, dim_type dim)
        void sign_normalize()
        void strong_normalize()
        cppbool check_strong_normalized()

    # Operators for constraint class
    Con operator=(Con &c)
    Con operator<(Linear_Expr e1, const Linear_Expr& e2)
    Con operator<(Var v1, Var v2)
    Con operator<(Linear_Expr e, const FLINT_Integer& n)
    Con operator<(FLINT_Integer n, Linear_Expr e)
    Con operator>(Linear_Expr e1, const Linear_Expr& e2)
    Con operator>(Var v1, Var v2)
    Con operator>(Linear_Expr e, FLINT_Integer n)
    Con operator>(FLINT_Integer n, Linear_Expr e)
    Con operator==(Linear_Expr e1, const Linear_Expr& e2)
    Con operator==(Var v1, Var v2)          
    Con operator==(Linear_Expr e, FLINT_Integer n)
    Con operator==(FLINT_Integer n, Linear_Expr e)
    Con operator<=(Linear_Expr e1, const Linear_Expr& e2)
    Con operator<=(Var v1, Var v2)
    Con operator<=(Linear_Expr e, FLINT_Integer n)
    Con operator<=(FLINT_Integer n, Linear_Expr e)
    Con operator>=(Linear_Expr e1, const Linear_Expr& e2)
    Con operator>=(Var v1, Var v2)
    Con operator>=(Linear_Expr e, FLINT_Integer n)
    Con operator>=(FLINT_Integer n, Linear_Expr e)
    Con operator<(Affine_Expr e1, const Affine_Expr& e2)
    Con operator<(Affine_Expr e, const FLINT_Integer& n)
    Con operator<(Affine_Expr e, Var v)
    Con operator<(const FLINT_Integer& n, Affine_Expr e)
    Con operator<(Var v, Affine_Expr e)
    Con operator>(Affine_Expr e1, const Affine_Expr& e2)
    Con operator>(Affine_Expr e, const FLINT_Integer& n)
    Con operator>(Affine_Expr e, Var v)
    Con operator>(const FLINT_Integer& n, Affine_Expr e)
    Con operator>(Var v, Affine_Expr e)
    Con operator==(Affine_Expr e1, const Affine_Expr& e2)
    Con operator==(Affine_Expr e, const FLINT_Integer& n)
    Con operator==(Affine_Expr e, Var v)
    Con operator==(const FLINT_Integer& n, Affine_Expr e)
    Con operator==(Var v, Affine_Expr e)
    Con operator<=(Affine_Expr e1, const Affine_Expr& e2)
    Con operator<=(Affine_Expr e, const FLINT_Integer& n)
    Con operator<=(Affine_Expr e, Var v)    
    Con operator<=(const FLINT_Integer& n, Affine_Expr e)
    Con operator<=(Var v, Affine_Expr e)
    Con operator>=(Affine_Expr e1, const Affine_Expr& e2)
    Con operator>=(Affine_Expr e, const FLINT_Integer& n)
    Con operator>=(Affine_Expr e, Var v)
    Con operator>=(const FLINT_Integer& n, Affine_Expr e)
    Con operator>=(Var v, Affine_Expr e)

    # TODO: Implement the below functions
    # Con complement_con(const Con&c, Topol t)
    # std::pair<Con, Con> integral_complement_eq(const Con& c)
    # std::pair<Con, Con> integral_complement_cons(const Con& c)
    # cppbool is_integral_inconsistent(const Con& c)

    # "pplite/Gen.hh"

    cdef cppclass Gen
    cdef cppclass Gen:
        # int compare(const Gen& x, const Gen&y)
        enum GenType "Gen::Type":
            LINE
            RAY
            POINT
            CLOSURE_POINT
        struct Impl "Gen::Impl":
            Linear_Expr expr
            FLINT_Integer inhomo
            GenType type
            Impl()
            Impl(Linear_Expr, e, FLINT_Integer d, GenType t)
        Gen()
        # Gen(const Gen& g)
        Gen(Gen& g)
        Gen(GenType t, Linear_Expr e, FLINT_Integer d)
        Impl& impl()
        GenType type()
        void set_type(GenType t)
        cppbool is_line()   
        cppbool is_ray()
        cppbool is_point()
        cppbool is_closure_point()
        cppbool is_line_or_ray()
        cppbool is_point_or_closure_point()
        Linear_Expr linear_expr()
        FLINT_Integer coeff(Var v)
        FLINT_Integer divisor()
        dim_type space_dim()
        void set_space_dim(dim_type dim)
        cppbool is_equal_to(const Gen& y)
        cppbool check_inv()
        void m_swap(Gen& y)
        # void permute_space_dims_cycle(const Dims& cycle, dim_type d)
        void shift_space_dims(dim_type start, dim_type n) 
        cppbool is_line_or_equality()
        void sign_normalize()
        void strong_normalize()
        cppbool check_strong_normalized()

    # Gen line(Linear_Expr e)
    # Gen ray(Linear_Expr e)
    # Gen point(Linear_Expr e, FLINT_Integer d) # double check
    # Gen closure_point(Linear_Expr e, FLINT_Integer d)
    # cppbool operator==(const Gen& x, const Gen& y)
    # cppbool operator!=(const Gen& x, const Gen& y)
    # void erase_space_dims(Gens& gs, Iter first, Iter last)
    # Index_Set invalid_lines(const Gens& gs)
    # Index_Set invalid_rays(const Gens& gs)
    # void rase_higher_dims(Gens& gs, dim_type d)
    # Gen materialize(const Indices& is, const Gens& gs)

    # "pplite/Itv.hh"
    cdef struct Itv:
        # enum Kind "Itv::Kind"
        # # enum Kind "Itv::Kind":
        # #     UNIVERSE
        # #     L_BOUNDED
        # #     U_BOUNDED
        # #     LU_BOUNDED
        # #     EMPTY 
        # Kind kind
        FLINT_Rational lb
        FLINT_Rational ub
        Itv(Spec_Elem s) #
        Itv& empty()
        Itv& universe()
        Itv& zero()
        cppbool check_inv()
        cppbool is_empty()
        cppbool is_universe()
        cppbool has_lb()
        cppbool has_ub()
        cppbool inf_lb()
        cppbool inf_ub()
        cppbool is_bounded()
        cppbool is_singleton()
        cppbool is_zero()
        cppbool is_disjoint_from(Itv& y)
        cppbool intersects(Itv& y)
        size_t hash()
        FLINT_Rational length()
        dim_type num_min_cons()
        dim_type num_rays()
        cppbool contains(Itv& y)
        cppbool contains(FLINT_Integer& num, FLINT_Integer& den)
        cppbool operator==(Itv& y)
        void set_empty()
        void set_universe()
        void set_zero()
        void set_lb(FLINT_Rational value)
        void set_ub(FLINT_Rational value)
        void set_singleton(FLINT_Rational value)
        void unset_lb()
        void unset_ub()
        cppbool glb_assign(const Itv& y)
        void lub_assign(const Itv& y)
        void widen_assign(const Itv& y)
        cppbool refine_as_integral()
        void complement_assign()
        void add_assign(const Itv& y)
        void mul_assign(const FLINT_Rational& r)

    Itv itv_from_con_inhomo(const Con& c)
    Itv itv_from_itv_con(const Con& c)
    Itv split_itv(Itv& itv, const Con&c, cppbool integral)
    Con get_lb_con(Var var, const Itv& itv)
    Con get_ub_con(Var var, const Itv& itv)
    Con get_eq_con(Var var, const Itv& itv)


# "pplite/Bits.hh"
    
    # cppclass Bits
    # cppclass Bits:
    #     Bits


    # cdef cppbool KVI
    # value template
    # ctypedef cppbool KVI "BBox.hh::KVI"
    # cdef cppclass Box
    # cdef cppclass Box:
    #     Box(Box& y)
    # #     Box& operator=(Box& y)
    # #     # Box(const Box<not KVI>& y)
    #     Box()
    #     Box(dim_type sd, Spec_Elem se)
    #     Volume_Info compute_volume_info()
    #     void maybe_update_volume_info()
    #     dim_type num_rays()
    #     FLINT_Rational pseudo_volume()
    #     cppbool check_inv()
    #     dim_type space_dim()
    #     dim_type affine_dim()
    #     cppbool is_empty()
    #     cppbool is_universe()
    #     cppbool is_bounded()
    #     cppbool is_bounded_expr()
    #     Topol topology()
    #     cppbool is_topologially_closed()
    #     cppbool constrains(Var v)
    #     Itv get_bounds(Var v)
    #     cppbool is_included_in(Con& c)
    #     cppbool contains(Box& y)
    #     cppbool equals(Box& y)
    #     cppbool is_disjoint_from(Box& y)
    #     cppbool less(Box& y)
    #     dim_type num_min_cons()
    #     dim_type num_min_gens()
    #     # Gens_Info gens_info()
    #     size_t hash()
    #     void swap(Box& y)
    #     void set_empty()
    #     void set_origin()
    #     void affine_image(Var var, Linear_Expr expr, FLINT_Integer inhomo, FLINT_Integer den)
    #     void concatenate_assign(Box& y)
    #     void glb_assign(Box& y)
    #     void lub_assign(Box& y)
    #     void time_elapse_assign(Box& y)
    #     void unconstrain(Var var)
    #     # void unconstrain(Iter first, Iter last)      
    #     # void unconstrain(const Index_Set& vars)
    #     void widening_assign(const Box& y)
    #     void add_space_dims(dim_type d, cppbool project)
    #     # void permute_space_dims(const Dims& perm)
    #     # void map_space_dims(const Dims& pfunc)
    #     void remove_space_dim(Var var)
    #     # void remove_space_dims(const Index_Set& vars)
    #     void remove_higher_space_dims(dim_type new_dim)
    #     void expand_space_dim(Var var, dim_type d)
    #     # void fold_space_dims(const Index_Set& vars, Var dst)
    #     void add_gen(const Gen& g)
    #     # void add_gens(const Gens& gs)
    #     void refine_bounds(dim_type dim, const Itv& y_itv)
    #     void refine_as_integral(dim_type dim)
    #     void refine_as_integral(dim_type first, dim_type last)

    # "pplite/BBox.hh" 
    # Note: We are explicitly defining the two possible Box classes.

    cdef cppclass Box_t "pplite::Box<true>": # namespace trick
        ctypedef pair[dim_type, FLINT_Rational] Volume_Info "Volume_Info"
        Box_t(Box_f& y)
        Box_t(const Box_t&& y)
        Box_t(dim_type sd, Spec_Elem se)
        Volume_Info compute_volume_info()
        void maybe_update_volume_info()
        dim_type num_rays()
        FLINT_Rational pseudo_volume()
        cppbool check_inv()
        dim_type space_dim()
        dim_type affine_dim()
        cppbool is_empty()
        cppbool is_universe()
        cppbool is_bounded()
        cppbool is_bounded_expr()
        Topol topology()
        cppbool is_topologially_closed()
        cppbool constrains(Var v)
        Itv get_bounds(Var v)
        cppbool is_included_in(Con& c)
        cppbool contains(const Box_t& y)
        # cppbool contains(const Box_f& y)
        cppbool equals(Box_t& y)
        # cppbool equals(Box_f& y)
        cppbool is_disjoint_from(const Box_t& y)
        # cppbool is_disjoint_from(const Box_f& y)
        cppbool less(const Box_t& y)
        # cppbool less(const Box_f& y)
        dim_type num_min_cons()
        dim_type num_min_gens()
        # Gens_Info gens_info()
        size_t hash()
        void swap(Box_t& y)
        void swap(Box_f& y)
        void set_empty()
        void set_origin()
        void affine_image(Var var, Linear_Expr expr, FLINT_Integer inhomo, FLINT_Integer den)
        void concatenate_assign(Box_t& y)
        # void concatenate_assign(Box_f& y)
        void glb_assign(Box_t& y)
        # void glb_assign(Box_f& y)
        void lub_assign(Box_t& y)
        # void lub_assign(Box_f& y)
        void time_elapse_assign(Box_t& y)
        # void time_elapse_assign(Box_f& y)
        void unconstrain(Var var)
        # void unconstrain(Iter first, Iter last)      
        # void unconstrain(const Index_Set& vars)
        void widening_assign(const Box_t& y)
        void widening_assign(const Box_f& y)
        void add_space_dims(dim_type d, cppbool project)
        # void permute_space_dims(const Dims& perm)
        # void map_space_dims(const Dims& pfunc)
        void remove_space_dim(Var var)
        # void remove_space_dims(const Index_Set& vars)
        void remove_higher_space_dims(dim_type new_dim)
        void expand_space_dim(Var var, dim_type d)
        # void fold_space_dims(const Index_Set& vars, Var dst)
        void add_gen(const Gen& g)
        # void add_gens(const Gens& gs)
        void refine_bounds(dim_type dim, const Itv& y_itv)
        void refine_as_integral(dim_type dim)
        void refine_as_integral(dim_type first, dim_type last)


    cdef cppclass Box_f "pplite::Box<false>":
        ctypedef pair[dim_type, FLINT_Rational] Volume_Info "Volume_Info"
        Box_f(Box_t& y)
        # Box_f(const Box_f&& y)
        Box_f(dim_type sd, Spec_Elem se)
        Volume_Info compute_volume_info()
        void maybe_update_volume_info()
        dim_type num_rays()
        FLINT_Rational pseudo_volume()
        cppbool check_inv()
        dim_type space_dim()
        dim_type affine_dim()
        cppbool is_empty()
        cppbool is_universe()
        cppbool is_bounded()
        cppbool is_bounded_expr()
        Topol topology()
        cppbool is_topologially_closed()
        cppbool constrains(Var v)
        Itv get_bounds(Var v)
        cppbool is_included_in(Con& c)
        cppbool contains(Box_t& y) 
        cppbool contains(Box_f& y)
        cppbool equals(Box_t& y)
        cppbool equals(Box_f& y)
        cppbool is_disjoint_from(Box_t& y)
        cppbool is_disjoint_from(Box_f& y)
        cppbool less(Box_t& y)
        cppbool less(Box_f& y)
        dim_type num_min_cons()
        dim_type num_min_gens()
        # Gens_Info gens_info()
        size_t hash()
        void swap(Box_t& y)
        void swap(Box_f& y)
        void set_empty()
        void set_origin()
        void affine_image(Var var, Linear_Expr expr, FLINT_Integer inhomo, FLINT_Integer den)
        void concatenate_assign(Box_t& y)
        void concatenate_assign(Box_f& y)
        void glb_assign(Box_t& y)
        void glb_assign(Box_f& y)
        void lub_assign(Box_t& y)
        void lub_assign(Box_f& y)
        void time_elapse_assign(Box_t& y)
        void time_elapse_assign(Box_f& y)
        void unconstrain(Var var)
        # void unconstrain(Iter first, Iter last)      
        # void unconstrain(const Index_Set& vars)
        void widening_assign(const Box_t& y)
        void widening_assign(const Box_f& y)
        void add_space_dims(dim_type d, cppbool project)
        # void permute_space_dims(const Dims& perm)
        # void map_space_dims(const Dims& pfunc)
        void remove_space_dim(Var var)
        # void remove_space_dims(const Index_Set& vars)
        void remove_higher_space_dims(dim_type new_dim)
        void expand_space_dim(Var var, dim_type d)
        # void fold_space_dims(const Index_Set& vars, Var dst)
        void add_gen(const Gen& g)
        # void add_gens(const Gens& gs)
        void refine_bounds(dim_type dim, const Itv& y_itv)
        void refine_as_integral(dim_type dim)
        void refine_as_integral(dim_type first, dim_type last)