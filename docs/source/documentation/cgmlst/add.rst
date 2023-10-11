.. _cgmlst_add:

.. toctree::
    :glob:

===========================
Add strains to the database
===========================

Next, you need to iteratively add your strains to the database. You
can use a draft genome (we recommend using `Spades
<http://cab.spbu.ru/software/spades/>`_ for assembly).
You can also add a reference genome for comparison.


.. note::

   You need to add each strain one by one to the database. You can
   specify strain name using **-s** option.


Genome data
^^^^^^^^^^^

You can add strains using GENOME fasta sequence file.

.. code-block:: bash
   
   wgMLST add --help
   Usage: wgMLST add [OPTIONS] DATABASE GENOME
   
   Adds a strain GENOME to the wgMLST DATABASE.
   
   Options:
   -s, --strain TEXT     Name of the strain (default:genome name)
   -i, --identity FLOAT  Minimum identity to search gene (default=0.95)
   -c, --coverage FLOAT  Minimum coverage to search gene (default=0.9)


Reads data
^^^^^^^^^^

Alternatively, you can also add strains from raw reads direcly with
single or paired FASTQS(.gz) files.

.. code-block:: bash
   
   wgMLST add2 --help
   Usage: wgMLST add2 [OPTIONS] DATABASE [FASTQS]...
   
   Adds a strain from FASTQS(.gz) reads to the wgMLST DATABASE.
   
   Options:
   -s, --strain TEXT     Name of the strain (default:genome name).
   -i, --identity FLOAT  Minimum identity to search gene (default=0.95).
   -c, --coverage FLOAT  Minimum coverage to search gene (default=0.9).
   -r, --reads INTEGER   Minimum reads coverage to search a gene
                         (default=10).


.. note::

   The defaut identity and coverage treshold are set to 0.9 and can be
   modulated with **-i** and **-c** options.


.. warning::

   Carefully check that the allele calling has been performed
   correctly for each genome. Check the number of genes found for each
   strain using the :ref:`strain command <strain_check>`.
   
