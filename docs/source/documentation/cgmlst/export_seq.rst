.. _cgmlst_export_seq:

.. toctree::
    :glob:

================
Export sequences
================

You can access to allele sequences present in the database and specify
a list of genes to export with **-l** option.

.. note::

   The gene list can be obtained with the :ref:`gene <gene_check>` command.

.. _cgmlst_export_sequence:

Sequence
========

A simple export of the different sequences.

.. code-block:: bash

   wgMLST sequence -h
   Usage: wgMLST sequence [OPTIONS] DATABASE

   Extracts sequences from a wgMLST DATABASE.
   
   Options:
   -o, --output FILENAME  Output result in fasta format (default:stdout).
   -f, --file FILENAME    File containing list of coregenes to extract
                          (default:all coregenes).
   --reference            Returns reference sequence instead of strain alleles.

.. _cgmlst_export_msa:

MSA
===

A multialign fasta file with concatenated genes. The file can be used
directly for phylogenetic analysis using maximum likelihood or
Bayesian approaches.

.. code-block:: bash
   
   wgMLST msa -h
   Usage: wgMLST msa [OPTIONS] DATABASE
   
   Computes Multiple Sequence Alignment from a wgMLST DATABASE.
   
   Options:
   ...
   -r, --realign          Realigns genes with same length (Default:No).


.. warning::

   It is highly recommended to define a limited list of genes to be
   exported for the phylogenetic approach.
   

