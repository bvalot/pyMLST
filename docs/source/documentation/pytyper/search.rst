.. _pytyper_search:

.. toctree::
    :glob:

=======================================
Search other typing profile of a strain
=======================================

In complement to classical MSLT, you can search for other typing METHOD
using a draft genome or raw reads.


Typing method
^^^^^^^^^^^^^

fimH typing
   FimH typing is based on the allelic sequence of the fimH gene in the
   species *Escherichia coli*
   (`Dias et al, 2010 <https://journals.asm.org/doi/10.1128/jcm.01858-09>`_).
   Allelic sequence were download from `CGE <https://bitbucket.org/genomicepidemiology/fimtyper_db/>`_.

spa typing
   Spa typing is based on the repetitions polymorphism present on the protein
   A gene (spa) in the species *Staphylococcus aureus*
   (`Frénay et al, 1996 <https://link.springer.com/article/10.1007/BF01586186>`_).
   Repetitions and sequence types definition were download from
   `Ridom <https://spa.ridom.de/>`_.

Clermont typing
   Clermont phylogrouping is based on the presence/absence of 4 different genes
   in the species *Escherichia coli*
   (`Clermont et al, 2012 <https://enviromicro-journals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12019>`_).   


Genome data
^^^^^^^^^^^

You can search typing METHOD from GENOME fasta sequence files.

.. code-block:: bash

   pyTyper search --help
   Usage: pyTyper search [OPTIONS] {fim|spa|clmt} GENOMES...

   Searches strain type using specified METHOD for an assembly GENOME.

   fim: fimH typing for Escherichia coli
   spa: spa typing for Staphylococcus aureus
   clmt: Phylogouping using ClermontTyping method for Escherichia coli

   Options:
   -i, --identity FLOAT   Minimum identity to search gene.
   -c, --coverage FLOAT   Minimum coverage to search gene.
   -f, --fasta FILENAME   Writes fasta file with gene allele.
   -o, --output FILENAME  Writes search result to (default:stdout).


.. note::

   If new alleles are present or you want to have sequence target by the typing method in your strains,
   you can obtain their sequences with the **-f** option.

.. note::

   You can perform searches on multiple genomes simultaneously.
