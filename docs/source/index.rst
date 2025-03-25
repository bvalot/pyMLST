.. PyMLST documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

******
pyMLST
******

.. figure:: logo.png
   :align: center
   :height: 150px
   :alt: pyMLST

   python Mlst Local Search Tool

Purpose
=======


Bacterial typing is critical to unraveling the spread of pathogens.
For this purpose, data from next-generation sequencing are
now widely used, with core multilocus sequence typing (cgMLST) or
whole genome multilocus sequence typing (wgMLST) becoming the new
standard. These methods are an extension of the traditional MLST
method, which uses a short list of housekeeping genes. cgMLST and
wgMLST use a large set of genes corresponding to the core or whole
genome. Similar to MLST, each unique sequence corresponds to a
specific allele, and the combination of alleles determines the
sequence type (ST) of the strain.


We have developed pyMLST to perform this task. Unlike other tools, it
uses a local SQLite database to store allele sequences and MLST
profiles. This allows the collection of strains to be expanded
iteratively. The input can be (i) an assembler-generated draft
genome, (ii) the direct raw data, or (iii) other genomes stored in the
sequence database.


Documentation
=============

.. toctree::
   :maxdepth: 2
   :caption: Users:
   
   documentation/installation
   documentation/cgmlst
   documentation/clamlst
   documentation/pytyper

   
.. toctree::
   :maxdepth: 2
   :caption: Developers:

   development
   api


Citation
========

If you use pyMLST, please cite the following paper:

Bignenet A. et al., Introduction and benchmarking of pyMLST:
open-source software for assessing bacterial clonality using core
genome MLST. 2023 Microbials Genomics, 9(11), 1126.
doi: `10.1099/mgen.0.001126 <https://doi.org/10.1099/mgen.0.001126>`_
