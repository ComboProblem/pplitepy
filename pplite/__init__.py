r"""
Cython wrapper for the Parma Polyhedra Lite Library (PPLite).

The Parma Polyhedra Lite Library (PPLite) is a library for polyhedral computations over
the rationals. This interface tries to reproduce the C++ API as faithfully as possible
in Python.

AUTHORS:

- Acadia Larsen (2024): initial version.
"""

__version__ = "0.0.51"

from .linear_algebra import (
        Variable, Linear_Expression, Affine_Expression
        )

from .constraint import (
        Constraint
        )

from .generators import (
        PPliteGenerator, Point, Closure_point, Ray, Line
        )

from .intervals import (
        Interval
        )

from .bounding_box import (
        Bounding_Box_t, Bounding_Box_f
        )

from .polyhedron import (
        NNC_Polyhedron
        )