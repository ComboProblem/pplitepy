# # distutils: language = c++
# # distutils: libraries = gmp gmpxx pplite m flint

# cimport cython

# from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check

# from libcpp.vector cimport vector as cppvector

# import_gmpy2()

# cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
#     r""" Converts FLINT_Integer to python integer.
    
#     TESTS::

#         >>> cdef fmpz_t x
#         >>> fmpz_init(x)
#         >>> fmpz_set_ui(x, 7)
#         >>> cdef FLINT_Integer &w
#         >>> w = new FLINT_Integer(x)
#         >>> fmpz_clear(x)
#         >>> z = FLINT_Integer_to_Python(w)
#         >>> print(z)
#         7
    
#     """
#     cdef mpz_t new_int
#     mpz_init(new_int)
#     fmpz_get_mpz(new_int, integer.impl())
#     y = GMPy_MPZ_From_mpz(new_int)
#     mpz_clear(new_int)
#     return y
    



    
# def integer_to_python_integer(integer):
#     if isinstance(fmpz_t, integer):
#         cdef mpz
# ctypedef fused fmpz_or_mpz:
#     fmpz
#     mpz

# cpdef fmpz_or_mpz integer_to_python(fmpz_or_mpz integer):
#     if isinstance(integer, fmpz):
#         pass
#     return GMPy_MPZ_From_mpz(integer).get_mpz_t()

# def integer_to_python(obj):
#     if isinstance(obj, fmpz): #no 
#         fmpz(obj)
#     if isinstance(obj, mpz):
#         return GMPy_MPZ_From_mpz(obj).get_mpz_t()

# def Integer_vector_to_list(Integer &cpp_vector):
#     cdef int i
#     cdef int size = Integer.size()
#     cdef list Integer_List = []
#     for i in range(size):
#         Integer_List.append(Integer[i])
#     return Integer_List


# cdef  PPL_Coefficient_from_pyobject(c) except *:
#     cdef mpz coeff

#     if MPZ_Check(c):
#         coeff = <mpz> c
#     elif isinstance(c, (int, str)):
#         coeff = mpz(c)
#     else:
#         coeff = mpz(c)
#         if coeff != c or c != coeff:
#             raise TypeError('ppl coefficients must be integral')

#     return PPL_Coefficient(coeff.z)

# def test_fun():
# ## TEST BLOCK BUT NOT PROPPER
#     #fmpz_t x, y;
#     #fmpz_init(x);
#     #fmpz_init(y);
#     # fmpz_set_ui(x, 7);
#     # fmpz_mul(y, x, x);
#     # fmpz_print(x);
#     # flint_printf("^2 = ");
#     # fmpz_print(y);
#     # flint_printf("\n");
#     # fmpz_clear(x);
#     # fmpz_clear(y);
#     cdef fmpz_t x
#     fmpz_init(x)
#     fmpz_set_ui(x, 7)
#     fmpz_print(x)
# # print(isinstance(mpz_t, test_integer))
# # print(isinstance(fmpz_t, test_integer))
#     # flint_to_gmp(x)
#     cdef mpz_t new_int
#     mpz_init(new_int)
#     fmpz_get_mpz(new_int, x)
#     fmpz_clear(x)
#     y = GMPy_MPZ_From_mpz(new_int)
#     mpz_clear(new_int)
#     print(y)
    
    
#     # print("We've read some code")
# # print()
# # print(isinstance(mpz, changed_int))

# def integer_to_python():
#     pass


