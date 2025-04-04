.. _cgmlst_initialise:

.. toctree::
    :glob:

=====================
Initialise a database
=====================

The first step of a cg/wgMLST analysis is to initialise a database
by a list of genes with a reference sequence for each of them.

:cgMLST: A list of genes corresponding to the coregenome of a species.

:wgMLST: A list  of genes corresponding to the whole genome of a
		 species or a clone.

Import from cgmlst.org
======================

You can automatically import a cgMLST resource from `cgmlst.org
<https://www.cgmlst.org/>`_.

.. code-block:: bash
				
   wgMLST import -h
   Usage: wgMLST import [OPTIONS] DATABASE [SPECIES]...

   Creates a wgMLST DATABASE from an online resource.
   
   The research can be filtered by adding a SPECIES name.
   
   Options:
   -f, --force             Overwrite alrealdy existing DATABASE
   --prompt / --no-prompt  Do not prompt if multiple choices are found,
				           fail instead.


Create from external scheme
===========================

The cg/wgMLST database can be created using a **scheme** corresponding to
a list of different genes (a multi-fasta file containing gene
sequences in nucleotide format).

.. code::

   >ACICU_RS02500
   TTATTTCTTCACAACAGATGGTGCAATTGGGTCGGCAGTGATATAGCCAACTGCTGCTGC
   ...
   GTGGTTAGAAGCAGTGGTCAT
   >ACICU_RS11305
   CGCACCTAATGGAAGAAAAGGGATCCCCGTAAACCATTTTAAAATATCGCGACGTGTTGG
   ...
   TTTGGAATTGATGCAGAAATTAAATCTTAA
   >ACICU_RS08820
   ATGGCTTATCAAACTTTAGAACAGCTACAGCAGTCTAAAGCCAAGCTTCACGAAACTGTG
   ...
   TCGCAGTTACGTTAA

.. warning::

   At contrary to other cg/wgMLST tools, only one allele for each
   gene must be include on the scheme file.


You can get scheme for:

:cgMLST:

* Using a scheme from a scientific publication and not available on
  `cgmlst.org <https://www.cgmlst.org/>`_.
	 
* Using the annotation of the genes from the reference genome of
  the species. After adding your strains to the database, you can
  filter to core genome by removing genes absent from least 95% of
  the strains (see :ref:`validate <m_option_check>`)
	 
:wgMLST:

* Using gene annotations from a genome close to your strains

* Using pangenome results from analysis of your strains with
  e.g. `Roary <https://sanger-pathogens.github.io/Roary/>`_.
	
	   

.. code-block:: bash

   wgMLST create --help
   Usage: wgMLST create [OPTIONS] DATABASE COREGENE
   
   Creates a wgMLST DATABASE from a template COREGENE.
   
   Options:
   -f, --force        Overwrite alrealdy existing DATABASE
   -c, --concatenate  Automatically concatenates genes with duplicated sequences
   -r, --remove       Automatically removes genes with duplicated sequences
   -s, --species TEXT  Name of the species
   -V, --version TEXT  Version of the database


.. warning::

   If the same sequence is used more than once in your scheme, you can
   specify how to handle it using the **-c** or **-r** options.
   

