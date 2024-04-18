{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------
The easiest way to install `hydrant` is directly from the source code.
The source code is available via its GitHub repository's
`release`_ webpage. Additionally, stable versions are also accessible
Python's Package Index (`PyPI`_).

.. _release: https://github.com/kasra-keshavarz/hydrant/releases
.. _PyPI: https://pypi.org/project/hydrantpy/

Instructions for installing :ref:`from source <install.source>`,
:ref:`PyPI <install.pypi>`, or its
:ref:`development version <install.dev>` are provided in the following.

.. _install.version:

Python versions supported
-------------------------

Officially Python 3.10, 3.11 and 3.12.

Installing hydrant
------------------

.. _install.source:

Installing from source 
~~~~~~~~~~~~~~~~~~~~~~

The easiest way to install hydrant is to download a release from
its `release`_ webpage and install it using Python's native ``pip``
installer.


After downloading an official release, a typical installation may
proceed with:

.. code-block:: shell

   foo@bar:~ $ tar -xvf hydrant-|version|.tar.gz
   foo@bar:~ $ cd hydrant-|version|
   foo@bar:hydrant-|version| $ pip install .

After installation, hydrant could be directly imported via Python's
interpreter.

.. _install.pypi:

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

The hydrant package is also accessible from PPyPI via pip:

.. code-block:: shell
  :linenos:

  pip install hydrantpy

.. note::

  The package's name under the PyPI repositories is ``hydrantpy``

.. _install.source:

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

See the :ref:`contributing guide <contributing>` for complete instructions on building from the git source tree.
Further, see :ref:`creating a development environment <contributing_environment>` if you wish to create
a hydrant development environment.

Running the test suite
----------------------

Unit tests are under active development...


.. _install.dependencies:

Dependencies
------------

.. _install.required_dependencies:

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

pandas requires the following dependencies.

================================================================ ==========================
Package                                                          Minimum supported version
================================================================ ==========================
`GeoPandas <https://geopandas.org>`__                            0.14.3
`NetworkX <https://networkx.org>`__                              3.2.1
================================================================ ==========================

