from __future__ import absolute_import

from .pplite_decl cimport *

cdef class Variable:
    cdef Var *thisptr

cdef class Linear_Expression:
    cdef Linear_Expr *thisptr

cdef class Affine_Expression:
    cdef Affine_Expr *thisptr