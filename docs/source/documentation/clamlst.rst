.. _clamlst:

***********************
classical MLST analysis
***********************

A workflow analysis of classical MLST is performed using a series of
Python scripts described below.


.. toctree::
   :maxdepth: 1

   clamlst/initialise
   clamlst/search

All avalaible commands can be listed using help fonction:

.. code-block:: bash
   
   claMLST --help
   
   Usage: claMLST [OPTIONS] COMMAND [ARGS]...
   
   Classical MLST commands.
   
   Commands:
   create   Creates a classical MLST DATABASE from a SCHEME csv and ALLELES...
   import   Creates a claMLST DATABASE from an online resource.
   info     Output the informations about a classical MLST DATABASE
   remove   Removes ALLELE sequence from the GENE on a mlst DATABASE.
   search   Searches ST number for an assembly GENOME using an mlst DATABASE.
   search2  Searches ST number from FASTQS(.gz) raw reads using an mlst...
