.. _installation:

.. toctree::
    :glob:

============
Installation
============

This section provides instructions for installation and configuration of pyMLST.


Automatic Installation
======================

You can install pyMLST and their dependancy using `bioconda <https://anaconda.org/bioconda/pymlst>`_:

.. code-block:: bash

   conda install -c bioconda -c conda-forge pymlst

Manual Installation
===================

* From `pypi repository <https://pypi.org/project/PyMLST/>`_:

  .. code-block:: bash

	 pip install pymlst

* From `github source <https://github.com/bvalot/pyMLST/>`_:

  .. code-block:: bash

	 virtualenv venv
	 source venv/bin/activate
	 make install
	 make build


Dependancy
==========

PyMLST uses 3 external tools to run alignment:

* Mafft (>=7.307)
  
  .. code-block:: bash

	 sudo apt install mafft
	 
* Blat (v35). You need to compile source or obtain executable at:
  https://genome.ucsc.edu/FAQ/FAQblat.html
  
* kma (>=1.3) You need to compile source from:
  https://bitbucket.org/genomicepidemiology/kma/src/master/


Configuration
=============

Configure the executable locations (if they are not on the PATH) and log level :

.. code-block:: bash
				
   pyMLST configure --help
   Usage: pyMLST configure [OPTIONS]

   Configure executables paths and log level.

   Options:
   -b, --blat FILE   Blat executable absolute path.
   -k, --kma FILE    Kma executable absolute path.
   -m, --mafft FILE  Mafft executable absolute path.
   -l, --log [DEBUG|INFO|WARNING|ERROR]
                     Level of logging, default=INFO  
   -r, --reset       Reset the configuration.
   --help            Show this message and exit.



