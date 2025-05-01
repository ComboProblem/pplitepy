# distutils: extra_compile_args = -std=c++11

from __future__ import absolute_import

from .pplite_decl cimport *

cdef class NNC_Polyhedron(object):
    cdef Poly* thisptr

cdef class Polyhedron_Constraint_Rel(object):
    cdef Poly_Con_Rel* thisptr

cdef class Polyhedron_Generator_Rel(object):
    cdef Poly_Gen_Rel* thisptr



# cdef class Constraint_Proxy(object):
#     cdef Cons_Proxy c_proxy