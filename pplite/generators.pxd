# distutils: extra_compile_args = -std=c++11

from __future__ import absolute_import

from .pplite_decl cimport *

from .constraint cimport *

cdef class PPliteGenerator(object):
    cdef Gen *thisptr

cdef enum GenType:
    LINE
    RAY
    POINT
    CLOSURE_POINT