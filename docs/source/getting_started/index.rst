{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------
The easiest way to install hydrant is to directly install it 
from the relevant GitHub repository. The package is also accessible
from the PyPI with the name ``hydrantpy``. 

Instructions for installing :ref:`from source <install.source>`,
:ref:`PyPI <install.pypi>`, or a
:ref:`development version <install.dev>` are also provided.

.. _install.version:

Python versions support
-----------------------

Officially Python 3.9, 3.10, 3.11 and 3.12.

Installing hydrant
------------------

.. _install.github:

Installing from GitHub 
~~~~~~~~~~~~~~~~~~~~~~

The easiest way to get the most up-to-date version of hydrant
is using Python's native ``pip`` installer.

.. code-block:: shell

   pip install git+https://github.com/kasra-keshavarz/hydrant.git

The hydrant package could be imported by Python.

.. _install.pypi:

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

The hydrant package is also accessible from the PyPI via pip:

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

