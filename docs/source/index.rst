.. baseq documentation master file, created by
   sphinx-quickstart on Fri May 11 23:16:22 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BeiSeq
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Welcome to BeiSeq's documentation.

Install
--------
.. code-block:: sh

   pip install beiseq

Guidance For Developer
----------------------
.. toctree::
   :maxdepth: 2

   guidance/phil
   guidance/cmd
   guidance/docs

Develop Interfaces
------------------

.. toctree::
   :maxdepth: 2

   BAM
   CNV
   SNV
   RNA/index
   Drops/Drops

* :ref:`modindex`

Table
-----

=====  =====  ======
   Inputs     Output
------------  ------
  A      B    A or B
=====  =====  ======
False  False  False
True   False  True
False  True   True
True   True   True
=====  =====  ======

`Docs for this project <http://packages.python.org/an_example_pypi_project/>`