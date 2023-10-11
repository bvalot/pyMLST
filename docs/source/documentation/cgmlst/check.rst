.. _cgmlst_check:

.. toctree::
    :glob:

=============================
Check quality of the database
=============================

After :ref:`loading <cgmlst_add>` all your strains to the database,
you need to check allele calling quality before :ref:`export results
<cgmlst_export_res>`.

.. note::

   You can have information of current data in the database using
   **stats** command.

   .. code-block:: bash

	  wgMLST stats -h
	  Usage: wgMLST stats [OPTIONS] DATABASE

	  Extract stats from a wgMLST DATABASE.


.. _strain_check:

Validate strains
================

To search potential strain with problems like bad assembly or wrong
species, you can use the **strain** command with the **-c** option.


.. code-block:: bash
				
   wgMLST strain -h
   Usage: wgMLST strain [OPTIONS] DATABASE

   Extracts a list of strains from a wgMLST DATABASE.
   
   Options:
   -m, --mincover INTEGER  Minimun number of strain found to keep a gene
                           (default:0)
   -k, --keep              Keep only gene with different allele (omit missing).
   -d, --duplicate         Conserve duplicate gene (default remove).
   -V, --inverse           Keep only gene that do not meet the filter
                           of mincover or keep options.
   -c, --count             Count the number of gene present in the database for
                           each strains.
   -o, --output FILENAME   Export strain list to (default=stdout).

..  note::
	
	If some strains show low number of genes found in comparison to the
	other, you can remove it using :ref:`remove <remove_check>`
	command.

.. note::

   Similarly to :ref:`gene <gene_check>` command or :ref:`export <cgmlst_export_res>`, you can filter gene
   that you want to conserved for the search.

   By default, only duplicate genes are removed.

.. _gene_check:

Validate genes
==============

Similarly to strains, it could be interesting to saved genes list to
conserved for the rest of the analysis using **gene** command.

.. code-block:: bash
				
   wgMLST gene -h
   Usage: wgMLST gene [OPTIONS] DATABASE
   
   Extracts a list of genes from a wgMLST DATABASE.
   
   Options:
   -m, --mincover INTEGER  Minimun number of strain found to keep a gene
				           (default:0)
   -k, --keep              Keep only gene with different allele (omit missing).
   -d, --duplicate         Conserve duplicate gene (default remove).
   -V, --inverse           Keep only gene that do not meet the filter of
                           mincover or keep options.
   -o, --output FILENAME   Export GENE list to (default=stdout).

.. note::

   Gene list that pass your threshold can be used further for :ref:`export
   sequence <cgmlst_export_seq>`. 
   
.. _m_option_check:
   
.. warning::

   An important parameter are the **-m** option that defined the
   minimum number of strains found to keep a gene.

   If you are interesting by coregene, you can defined this number to
   correspond to **95%** of the strain in the database.
   (As example, if you have 100 strains in your database, you need to
   set this parameter to 95)
   
   
.. _remove_check:

Remove strains or genes
=======================

After checking the database, if some strains or genes need to be
removed, you can use the **remove** commands.

.. code-block:: bash
	  
   wgMLST remove -h
   Usage: wgMLST remove [OPTIONS] DATABASE [GENES_OR_STRAINS]...
	  
   Removes STRAINS or GENES from a wgMLST DATABASE.
   
   Options:
   --strains / --genes    Choose the item you wish to remove  [default: strains]
   -f, --file FILENAME    File list of genes or strains to removed on the wgMLST
				          database.
