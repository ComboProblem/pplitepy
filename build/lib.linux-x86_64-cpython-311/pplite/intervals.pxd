from __future__ import absolute_import

from .pplite_decl cimport *

cdef class Interval(object):
    cdef Itv interval

# cpdef object Interval():

# cdef enum Spec_Elem "":
#     EMPTY
#     UNIVERSE

# cdef enum Kind "Itv::Kind":
#     UNIVERSE
#     L_BOUNDED
#     U_BOUNDED
#     LU_BOUNDED
#     EMPTY
