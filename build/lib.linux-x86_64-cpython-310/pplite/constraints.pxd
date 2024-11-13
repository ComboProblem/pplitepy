from .linear_algebra cimport *

cdef _make_Constraint_from_richcmp(lhs_, rhs_, op)

cdef class Constraint(object):
	cdef Con *thisptr