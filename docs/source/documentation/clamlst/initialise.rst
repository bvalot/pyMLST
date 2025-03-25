.. _clamlst_initialise:

.. toctree::
   :glob:

==========================
Initialise a MLST database
==========================

A MLST database contains the different alleles for each gene of the
scheme and a table of association of the alleles to determined the
sequence type (ST).


Import from pubMLST
===================

You can automatically import a MLST resource from `pubmlst
<https://pubmlst.org/data/>`_ or `pasteur <https://bigsdb.pasteur.fr/>`_.

.. code-block:: bash

   claMLST import -h
   Usage: claMLST import [OPTIONS] DATABASE [SPECIES]...

   Creates a claMLST DATABASE from an online resource.

   The research can be filtered by adding a SPECIES name.

   Options:
   --prompt / --no-prompt  Do not prompt if multiple choices are found,
				           fail instead.
   -f, --force        	   Overwrites alrealdy existing DATABASE
   -m, --mlst TEXT         Specifies the desired MLST scheme name.
   -r, --repository        Choose the online repository to use [pubmlst|pasteur]



Create from other resource
==========================

Alternatively, you can create a database with the allele sequence and
MLST profile of your favorite species. To create a database, pyMLST
needs the gene name in the MLST profile header to match the name in
the fasta file. For example, the rpoB gene in the MLST profile header
must match the rpoB.fas file. You will also need to remove the
additional column corresponding to the clonal complex in the MLST
profile file, if present.

.. code-block:: bash
   
   claMLST create --help
   Usage: claMLST create [OPTIONS] DATABASE PROFILE ALLELES...

   Creates a classical MLST DATABASE from a txt PROFILE and fasta ALLELES files.

   Options:
   -f, --force        	  Overwrites alrealdy existing DATABASE
   -s, --species TEXT     Name of the species
   -V, --version TEXT     Version of the database


   
Scheme example
--------------

.. code::
   
   ST      cpn60   fusA    gltA    pyrG    recA    rplB    rpoB
   1       1       1       1       1       5       1       1
   2       2       2       2       2       2       2       2
   3       3       3       2       2       3       1       3
   ...
		  
Allele example
--------------

.. code::
   
   >cpn60_1
   ATGAACCCAATGGATTTAAAACGCGGTATCGACATTGCAGTAAAAACTGTAGTTGAAAAT
   ATCCGTTCTATTGCTAAACCAGCTGATGATTTCAAAGCAATTGAACAAGTAGGTTCAATC
   TCTGCTAACTCTGATACTACTGTTGGTAAACTTATTGCTCAAGCAATGGAAAAAGTAGGT
   AAAGAAGGCGTAATCACTGTAGAAGAAGGTTCTGGCTTCGAAGACGCATTAGACGTTGTA
   GAAGGTATGCAGTTTGACCGTGGTTATATCTCTCCGTACTTTGCAAACAAACAAGATACT
   TTAACTGCTGAACTTGAAAATCCGTTCATTCTTCTTGTTGATAAAAAAATCAGCAACATT
   CGTGAATTGATTTCTGTTTTAGAAGCAGTTGCTAAAACTGGTAAA
   >cpn60_2
   ATGAACCCAATGGATTTAAAACGCGGTATCGACATTGCAGTAAAAACTGTAGTTGAAAAT
   ATCCGTTCTATTGCTAAACCAGCTGATGATTTCAAAGCAATTGAACAAGTAGGTTCAATC
   TCTGCTAACTCTGATACTACTGTTGGTAAACTTATTGCTCAAGCAATGGAAAAAGTAGGT
   AAAGAAGGCGTAATCACTGTAGAAGAAGGCTCAGGCTTCGAAGACGCATTAGACGTTGTA
   GAAGGTATGCAGTTTGACCGTGGTTATATCTCTCCGTACTTTGCAAACAAACAAGATACT
   TTAACTGCTGAACTTGAAAATCCGTTCATCCTTCTTGTTGATAAAAAAATCAGCAACATT
   CGTGAATTGATTTCTGTTTTAGAAGCAGTTGCTAAAACTGGTAAA
   ...
