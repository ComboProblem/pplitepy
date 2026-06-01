# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check, mpq, MPQ_Check, GMPy_MPQ_From_mpz, MPZ
from libcpp.vector cimport vector as cppvector
from libc.stdlib cimport malloc, free
import sys


# It is assumed that the pplite environment is set up to use FLINT_integers.
# TODO:  Write a proper conversion module to handle the Integer class in PPLite so this works regardless of setup.

import_gmpy2()

####################################################
########## BEGIN ADAPTATION BLOCK ##################
####################################################

# Code here is adapted from python-flint.
# Relevant modules: flint.types.fmpz.pyx, flint.types.fmpz.pxd
# TODO: Once the C-API is exposed on python-flint, change this to python-flint imports.

cdef fmpz_get_intlong(fmpz_t x):
    """
    Convert fmpz_t to a Python int or long.
    """
    cdef slong size
    cdef ulong * words
    cdef int i
    cdef fmpz_struct xabs[1]
    if COEFF_IS_MPZ(x[0]):
        # Python from signed bytes is slow so we convert the absolute value.
        # we need 8 * sizeof(ulong) * size >= fmpz_bits(x)
        size = fmpz_bits(x) // (8 * sizeof(ulong)) + 1
        words = <ulong *>malloc(size * sizeof(ulong))
        if fmpz_sgn(x) == -1:
            fmpz_init(xabs)
            fmpz_abs(xabs, x)
            fmpz_get_ui_array(words, size, xabs)
            fmpz_clear(xabs)
        else:
            fmpz_get_ui_array(words, size, x)
        if is_big_endian:
            for i in range(size):
                words[i] = ulong_from_little_endian(<unsigned char *>(words + i))
        v = int.from_bytes((<char *>words)[:size * sizeof(ulong)], "little")
        if fmpz_sgn(x) == -1:
            v = -v
        free(words)
        return v
    else:
        return <slong>x[0]

####################################################
############ END ADAPTATION BLOCK ##################
####################################################

cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
    """ Converts PPlite::FLINT_Integer to python int.
    """
    cdef fmpz_t pplite_integer = integer.impl()
    return fmpz_get_intlong(pplite_integer)

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer):
    r"""Converts python int to a PPLite::FLINT_Integer."""
    cdef fmpz_t x
    fmpz_init(x)
    fmpz_set_python(x, integer)
    return FLINT_Integer(x)

# TODO: Right now rationals are not used in the code exposed. When python-flint exposes the C-API

cdef FLINT_Rational_to_Python(FLINT_Rational& rational):
    """Converts the Flint_Rational c++ class to a python object.

    INPUT:

    - rational: FLINT_Rational (c++)

    OUTPUT:

    - mpq

    """
    cdef mpz_t a
    cdef mpz_t b
    mpz_init(a)
    mpz_init(b)
    fmpq_get_mpz_frac(a , b, rational.impl())
    frac = GMPy_MPQ_From_mpz(a, b)
    mpz_clear(a)
    mpz_clear(b)
    return frac

cdef FLINT_Rational Python_float_to_FLINT_Rational(rational):
    """ Converts python float or fraction to a FLINT_Rational (c++).

    INPUT:

    - rational: object with method .as_integer_ratio()

    OUTPUT:

    FLINT_Rational (c++) 
    """
    cdef FLINT_Integer num
    cdef FLINT_Integer den
    try:
        numerator, denominator = rational.as_integer_ratio()
    except ValueError:
        raise ValueError("Rational Conversion Failed.")
    num = Python_int_to_FLINT_Integer(numerator)
    dem = Python_int_to_FLINT_Integer(denominator)
    return FLINT_Rational(num, dem)

def FLINT_Integer_Conversion_Check(possible_integer):
    """
    Checks a python object is convertible to a FLINT_Integer.

    Input: Object

    Output: Bool
    """
    if isinstance(possible_integer, (int, str)):
        return True
    if MPZ_Check(possible_integer):
        return True
    return False
