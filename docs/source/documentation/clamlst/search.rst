.. _clamlst_search:

.. toctree::
    :glob:

===============================
Search MLST profile of a strain
===============================

Similarly to wgMLST analysis, you need a draft genome or raw reads
data to find the MLST profile. 

.. note::

   You can perform MLST searches on multiple genomes or raw reads
   simultaneously.
   
  

Genome data
^^^^^^^^^^^

You can search ST from GENOME fasta sequence files.

.. code-block:: bash
   
   claMLST search --help
   Usage: claMLST search [OPTIONS] DATABASE GENOMES
   
   Searches ST number for assembly GENOMES using an mlst DATABASE
   
   Options:
   -i, --identity FLOAT   Minimum identity to search gene (default=0.9)
   -c, --coverage FLOAT   Minimum coverage to search gene (default=0.9)
   -f, --fasta FILENAME   Writes fasta file with gene allele
   -o, --output FILENAME  Writes ST search result to (default:stdout)


Reads data
^^^^^^^^^^

Alternatively, you can search ST directly from raw reads with single
or paired FASTQS(.gz) files.

.. code-block:: bash
   
   claMLST search2 --help
   Usage: claMLST search2 [OPTIONS] DATABASE [FASTQS]...
   
   Searches ST number from FASTQS(.gz) raw reads using an mlst DATABASE.
   
   Options:
   -i, --identity FLOAT   Minimum identity to search gene (default=0.9).
   -c, --coverage FLOAT   Minimum coverage to search gene (default=0.95).
   -r, --reads INTEGER    Minimum reads coverage to search gene (default=10).
   --paired / --single    Defines type of fastqs files.
   -f, --fasta FILENAME   Writes fasta file with gene allele.
   -o, --output FILENAME  Writes ST search result to (default:stdout).


.. note::

   The default identity and coverage thresholds are set to 0.9 and can
   be modulated using the **-i** and **-c** options.

.. note::

   If new alleles are present, you can obtain their sequences with
   the **-f** option.

