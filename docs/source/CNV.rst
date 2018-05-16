.. _CNV:

CNV Analysis
============

CNV analysis is for: XXXX

Config
--------------
.. code-block:: sh

    [CNV]
    bowtie2 = /mnt/gpfs/Database/softs/anaconda2/bin/bowtie2
    samtools = /mnt/gpfs/Database/softs/anaconda2/bin/samtools

    [CNV_ref_hg19]
    bowtie2_index = /mnt/gpfs/Database/ref/hg19/hg19
    dynamic_bin = /mnt/gpfs/Users/zhangxiannian/basematic/cnv/hg19.dynabin.txt

.. code-block:: sh

    mkdir myproject
    cd myproject
    python3 -m venv venv

Dependencies
------------

These distributions will be installed automatically when installing Flask.

* `Werkzeug`_ implements WSGI, the standard Python interface between
  applications and servers.
* `Jinja`_ is a template language that renders the pages your application
  serves.
* `MarkupSafe`_ comes with Jinja. It escapes untrusted input when rendering
  templates to avoid injection attacks.
* `ItsDangerous`_ securely signs data to ensure its integrity. This is used
  to protect Flask's session cookie.
* `Click`_ is a framework for writing command line applications. It provides
  the ``flask`` command and allows adding custom management commands.

.. _Werkzeug: http://werkzeug.pocoo.org/
.. _Jinja: http://jinja.pocoo.org/
.. _MarkupSafe: https://pypi.org/project/MarkupSafe/
.. _ItsDangerous: https://pythonhosted.org/itsdangerous/
.. _Click: http://click.pocoo.org/

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

These distributions will not be installed automatically. Flask will detect and
use them if you install them.

* `Blinker`_ provides support for :ref:`signals`.
* `SimpleJSON`_ is a fast JSON implementation that is compatible with
  Python's ``json`` module. It is preferred for JSON operations if it is
  installed.
* `python-dotenv`_ enables support for :ref:`dotenv` when running ``flask``
  commands.
* `Watchdog`_ provides a faster, more efficient reloader for the development
  server.

.. _Blinker: https://pythonhosted.org/blinker/
.. _SimpleJSON: https://simplejson.readthedocs.io/
.. _python-dotenv: https://github.com/theskumar/python-dotenv#readme
.. _watchdog: https://pythonhosted.org/watchdog/

Virtual environments
--------------------

Use a virtual environment to manage the dependencies for your project, both in
development and in production.

What problem does a virtual environment solve? The more Python projects you
have, the more likely it is that you need to work with different versions of
Python libraries, or even Python itself. Newer versions of libraries for one
project can break compatibility in another project.

Virtual environments are independent groups of Python libraries, one for each
project. Packages installed for one project will not affect other projects or
the operating system's packages.

Python 3 comes bundled with the :mod:`venv` module to create virtual
environments. If you're using a modern version of Python, you can continue on
to the next section.

If you're using Python 2, see :ref:`install-install-virtualenv` first.

.. _install-create-env:

Create an environment
~~~~~~~~~~~~~~~~~~~~~

Create a project folder and a :file:`venv` folder within:

.. code-block:: sh

    mkdir myproject
    cd myproject
    python3 -m venv venv

On Windows:

.. code-block:: bat

    py -3 -m venv venv

If you needed to install virtualenv because you are on an older version of
Python, use the following command instead:

.. code-block:: sh

    virtualenv venv

On Windows:

.. code-block:: bat

    \Python27\Scripts\virtualenv.exe venv

.. _install-activate-env:

Activate the environment
~~~~~~~~~~~~~~~~~~~~~~~~

Before you work on your project, activate the corresponding environment:

.. code-block:: sh

    . venv/bin/activate

On Windows:

.. code-block:: bat

    venv\Scripts\activate

Your shell prompt will change to show the name of the activated environment.