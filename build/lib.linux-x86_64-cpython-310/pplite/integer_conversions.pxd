from __future__ import absolute_import

from .pplite_decl cimport *
####################################################

cdef FLINT_Integer_to_Python(FLINT_Integer& integer)

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer)

cdef FLINT_Rational_to_Python(FLINT_Rational& rational)

cdef FLINT_Rational Python_float_to_FLINT_Rational(rational)