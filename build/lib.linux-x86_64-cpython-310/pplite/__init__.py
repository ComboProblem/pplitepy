r"""
Cython wrapper for the Parma Polyhedra Lite Library (PPLite)

"""

__version__ = "0.0.27"

from .linear_algebra import (
        Variable, Linear_Expression, Affine_Expression
        )

from .constraint import (
        Constraint
        )

from .generators import (
        PPliteGenerator
        )

from .intervals import (
        Interval
        )

from .bounding_box import (
        Bounding_Box_t, Bounding_Box_f
        )