from libcpp cimport bool as cppbool
from libcpp.vector cimport vector as cppvector
from gmpy2 cimport mpz 

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
cdef extern from "pplite/pplite.hh" namespace "pplite":
# cdef extern from "pplite/FLINT_Integer.hh" namespace "pplite":
    cdef cppclass FLINT_Integer
    cdef cppclass FLINT_Integer:
        FLINT_Integer()
        const fmpz_t& impl() 
        FLINT_Integer(const fmpz_t z)
        FLINT_Integer(signed int si)
        FLINT_Integer(const mpz_t z)
        
# cdef extern from "pplite/GMP_Integer.hh" namespace "pplite":
#     cdef cppclass GMP_Integer
#     cdef cppclass GMP_Integer:
#         GMP_Integer()
#         const mpz_t& impl() 
#         # GMP_Integer(const fmpz_t z)
#         GMP_Integer(signed int si)
#         GMP_Integer(const mpz_t z)



# cdef extern from "pplite/globals.hh" namespace "pplite":
    ctypedef size_t dim_type

# cdef extern from "pplite/Var.hh" namespace "pplite":
    cdef cppclass Var
    cdef cppclass Var:
#       Var()
        Var(dim_type i)
        dim_type id()
        dim_type space_dim()
    ctypedef struct Vars_Set

    
# cdef extern from "pplite/Linear_Expr.hh" namespace "pplite":
    ctypedef cppvector[FLINT_Integer] Impl
    cdef cppclass Linear_Expr
    cdef cppclass Linear_Expr:
        Linear_Expr()
        Linear_Expr(Var v)
        Linear_Expr(dim_type dim)
        Linear_Expr(Linear_Expr &e)
        Linear_Expr(const Linear_Expr &e, dim_type dim)
        Linear_Expr operator+(Linear_Expr &e)
        Linear_Expr operator-(Linear_Expr &e)
        Linear_Expr operator*(FLINT_Integer &c)
        dim_type id()  #methods? #linear_expr.hh lines 39-112
        dim_type space_dim()
        void set_space_dim(dim_type dim)
        void swap_space_dims(dim_type i, dim_type j)
        void shift_space_dims(dim_type start, dim_type n)
        void shift_space_dims(Var start, dim_type n)
        Impl impl()
        cppbool is_equal_to(Linear_Expr& y)

    