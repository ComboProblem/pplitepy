# distutils: extra_compile_args = -std=c++11

from __future__ import absolute_import

from cpython.long cimport PyLong_Check
from .pplite_decl cimport *
import sys

####################################################
########## BEGIN ADAPTATION BLOCK ##################
####################################################

# Code here is adapted from python-flint.
# Relevant modules: flint.types.fmpz.pyx, flint.types.fmpz.pxd
# TODO: Once the C-API is exposed on python-flint, change this to python-flint imports.

cdef int is_big_endian = int(sys.byteorder == "big")

cdef inline ulong ulong_from_little_endian(unsigned char *ptr):
    # Read a ulong from little-endian bytes
    cdef ulong w = 0
    for i in range(sizeof(ulong) // 8):
        w = (w << 8) | ptr[i]
    return w

cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef int overflow
    cdef slong longval
    cdef slong size
    cdef bytes b
    cdef ulong w
    cdef ulong *words
    cdef int i

    longval = pylong_as_slong(<PyObject*>obj, &overflow)
    if overflow:
        # make sure the sign bit fits
        # we need 8 * sizeof(ulong) * size > obj.bit_length()
        size = obj.bit_length() // (8 * sizeof(ulong)) + 1
        b = obj.to_bytes(sizeof(ulong) * size, "little", signed=True)
        # b is a local Python object, we access the internal pointer
        words = <ulong*>(<char *>b)
        if is_big_endian:
            for i in range(size):
                words[i] = ulong_from_little_endian(<unsigned char *>(words + i))
        fmpz_set_signed_ui_array(x, words, size)
    else:
        fmpz_set_si(x, longval)

cdef inline int fmpz_set_python(fmpz_t x, obj):
    if PyLong_Check(obj):
        fmpz_set_pylong(x, obj)
        return 1
    return 0

####################################################
############ END ADAPTATION BLOCK ##################
####################################################

cdef FLINT_Integer_to_Python(FLINT_Integer& integer)

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer)

cdef FLINT_Rational_to_Python(FLINT_Rational& rational)

cdef FLINT_Rational Python_float_to_FLINT_Rational(rational)
