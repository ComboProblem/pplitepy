from __future__ import absolute_import

from .pplite_decl cimport *

cdef class Variable:
    cdef Var *thisptr

cdef class Linear_Expression:
    cdef Linear_Expr *thisptr
    
# cdef PPL_Coefficient PPL_Coefficient_from_pyobject(c) except *
