PPlite Python Wrapper
=====================

This Python package provides a wrapper to the C++ package `pplite <https://github.com/ezaffanella/PPLite/>`__.

How it works
------------

The names of objects and methods are a spelled out version of the objects and methods in the library:

.. code:: python

	>>> import pplite
	>>> x = pplite.Variable(0)
	>>> y = pplite.Variable(1)
	>>> z = pplite.Variable(2)
	>>> cons = [x >= 0, y >= 0, z > 0, x + y + z == 1]
	>>> poly = pplite.NNC_Polyhedron(spec_elem = "universe", dim_type = 3, topology = "nnc")
	>>> poly.add_constraints(cons)
	>>> poly
	x0+x1+x2-1==0, x0>=0, x1>=0, -x0-x1+1>0
	>>> poly.generators()
	[p(x2), c(x0), c(x1)]
	
The available objects and functions from `pplite` Python module are:

- `Variable`, `Linear_Expression`, `Affine_Expression` (defined in `pplite.linear_algebra`)

- `Constraint` (defined in `pplite.constraint`)

- `Interval` (defined in `pplie.intervals`)

- `Bounding_Box_t`, `Bounding_Box_f` (defined in `pplite.bounding_box`)

- `PPliteGenerator`, `Point`, `Closure_point`, `Line`, `Ray` (defined in `pplite.generators`)
	
- `NNC_Polyhedron` (defined in `pplite.polyhedron`)

Installation
------------

The project is available at `Python Package Index <https://pypi.org/project/pplitepy/>`_ and
can be installed with pip. It is reccomened to create a virtual enviorment
(suggested name `venv-pplite`) when installing the package.::

    $ python3 -m venv venv-pplite
    $ . ./venv-pplite/bin/activate
    (venv-pplite) $ pip install pplitepy
    
The package should be accessed from the virtual enviroment::

    (venv-pplite) $ python3
    >>> import pplite

Using from Cython
-----------------

All Python classes from pplpy are extension types and can be used with Cython. Most
extension type carries an attribute `thisptr` that holds a pointer to
the corresponding C++ object from ppl.

Source
------

You can find the latest version of the source code on github:
https://github.com/ComboProblem/pplitepy

Documentation
-------------

TBD.

License
-------

pplitepy is distributed under the terms of the GNU General Public License (GPL)
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version. See http://www.gnu.org/licenses/.

Requirements
------------

- `pplite <https://github.com/ezaffanella/PPLite>`__

- `Cython <http://cython.org>`_ (tested with 3.0)

- `cysignals <https://pypi.org/project/cysignals/>`_

- `gmpy2 <https://pypi.org/project/gmpy2/>`_

- `flint <https://flintlib.org/>`_

On Debian/Ubuntu systems the dependencies can be installed with::

    $ sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev libppl-dev cython3 python3-gmpy2 python3-cysignals-pari flint
