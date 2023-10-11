.. _make:

.. _using-the-makefile:

Using the `Makefile`
====================

This project includes a `Makefile <https://www.gnu.org/software/make/>`_
that you can use to perform common tasks such as running tests and building
documentation.

Targets
-------

This section contains a brief description of the targets defined in the
``Makefile``.

``clean``
^^^^^^^^^

Removes generated packages, documentation, temporary files, *etc*.

.. _make_lint:

``lint``
^^^^^^^^

Runs `pylint <https://www.pylint.org/>`_ against the project files.

.. _make_test:

``test``
^^^^^^^^

Runs the unit tests.

``quicktest``
^^^^^^^^^^^^^

Runs the unit tests without performing pre-test validations (like
:ref:`linting <make_lint>`).

.. _make_docs:

``docs``
^^^^^^^^

Builds the documentation for production.

.. note::

    You can also build the documents directly, bypassing validations like
    :ref:`linting <make_lint>` and :ref:`testing <make_test>` using
    `Sphinx Makefile <https://github.com/mapnik/sphinx-docs/blob/master/Makefile>`_
    directly.

    .. code-block:: bash

        cd docs
        make clean && make html
        make latexpdf

.. _make_answers:

``answers``
^^^^^^^^^^^

Performs a quick build of the documentation and open it in your browser.

``package``
^^^^^^^^^^^

Builds the package for publishing.

.. _make-publish:

``publish``
^^^^^^^^^^^

Publishes the package to your repository.

``build``
^^^^^^^^^

Installs the current project locally so that you may run the command-line application.

``venv``
^^^^^^^^

Creates a virtual environment.

``install``
^^^^^^^^^^^

Installs (or updates) project dependencies.

``install_docs``
^^^^^^^^^^^^^^^^

Installs (or updates) documentation dependencies.

``licenses``
^^^^^^^^^^^^

Generates a report of the projects dependencies and respective licenses.

.. note::

    If project dependencies change, please update this documentation.
