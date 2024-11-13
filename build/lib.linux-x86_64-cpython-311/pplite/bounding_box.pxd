from .pplite_decl cimport *

cdef class Bounding_Box_t:
    cdef Box_t *thisptr

cdef class Bounding_Box_f:
    cdef Box_f *thisptr

# cpdef Create_Bounding_Box(cppbool keep_volume_info, dim_type sd=*, Spec_Elem se=*)