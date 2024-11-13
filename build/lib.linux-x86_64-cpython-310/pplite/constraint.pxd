from .linear_algebra cimport *

cdef _wrap_Constraint(Con constraint)

cdef _make_Constraint_from_richcmp(lhs_, rhs_, op)

cdef class Constraint(object):
	cdef Con *thisptr

cdef enum ConType:
	EQUALITY
	NONSTRICT_INEQUALITY
	STRICT_INEQUALITY