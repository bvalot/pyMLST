.. _cgmlst_export_res:

.. toctree::
    :glob:

==============
Export results
==============

When the database is complete and :ref:`validated <cgmlst_check>`, you
can export results for futher analysis.

.. _cgmlst_export_distance:

Distance
========

A matrix of cgMLST distances can be computed from the database and
defined as the number of different alleles between each pair of two
strains, omitting the missing data.

.. code::
   
   #Strain  33_PA   34_PA   35_PA   61_PA   84_PA   98_PA 
   33_PA     0       39      37      25      20      23   
   34_PA     39      0       5       33      35      39   
   35_PA     37      5       0       31      33      37   
   61_PA     25      33      31      0       22      27   
   84_PA     20      35      33      22      0       21   
   98_PA     23      39      37      27      21      0  

.. code-block:: bash
				
   wgMLST distance --help
   Usage: wgMLST distance [OPTIONS] DATABASE
   
   Extracts a distance matrix from a wgMLST DATABASE.
   Options:
   -m, --mincover INTEGER  Minimun number of strains found to retain a gene
                           (default:0)
   -k, --keep              Keeps only gene with different alleles (omit
                           missing).
   -d, --duplicate         Keeps duplicate genes (default remove).
   -V, --inverse           Keeps only gene that do not match the filter of
                           mincover or keep options.
   
.. warning::

   To have correct distance calculation, you need to filter genes with
   low frequency observations. See :ref:`validate <m_option_check>`  to
   have more informations on **-m** option.

.. _cgmlst_export_mlst:
  
MLST
====

The MLST profiles can be also exported. The number indicated the
allele *id* in the database. An formatted version compatible with grapetree
can be defined.

.. code::
   
   #GeneId 33_PA   34_PA   35_PA   61_PA   84_PA   98_PA
   PA0120  3918    3918    3918    3918    3918    3918 
   PA0527  3963    3963    3963    3963    3963    3963 
   PA0691  3954    3954    3954    8945    3954    3954
   PA0935  3910    3910    3910    3910    3910    3910
   ...


.. code-block:: bash
				
   wgMLST mlst --help
   Usage: wgMLST mlst [OPTIONS] DATABASE

   Extracts an MLST table from a wgMLST DATABASE.
   Options:
   ...
   -f, --form [default|grapetree]  Specify format of output

.. note::

   Similarly to :ref:`distance <cgmlst_export_distance>`, the gene export on this mlst table can be
   defined with -m, -k, and -d options.

