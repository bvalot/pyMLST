[![PyPI version](https://badge.fury.io/py/PyMLST.svg)](https://pypi.org/project/PyMLST/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pymlst/README.html)
[![Documentation Status](https://readthedocs.org/projects/pymlst/badge/?version=latest)](https://pymlst.readthedocs.io/en/latest/?badge=latest)

# pyMLST
![pyMLST](docs/source/logo.png "A Python Mlst Local Search Tool")

A Python Mlst Local Search Tool.

## Purpose
Bacterial typing is critical to unraveling the spread of pathogens.
For this purpose, data from next-generation sequencing are now widely used, with core multilocus sequence typing (cgMLST) or whole genome multilocus sequence typing (wgMLST) becoming the new standard.
These methods are an extension of the traditional MLST method, which uses a short list of housekeeping genes.
cgMLST and wgMLST use a large set of genes corresponding to the core or whole genome.
Similar to MLST, each unique sequence corresponds to a specific allele, and the combination of alleles determines the sequence type (ST) of the strain.

We have developed pyMLST to perform this task. Unlike other tools, it uses a local SQLite database to store allele sequences and MLST profiles.
This allows the collection of strains to be expanded iteratively. The input can be (i) an assembler-generated draft genome, (ii) the direct raw data, or (iii) other genomes stored in the sequence database.

## New version
V2.1:

- Use raw reads (FASTQ) directly with the kma integration (search2 and add2 command)

V2.0:

- An automatic import database mechanism to initiated cgMLST and MLST databases.
- A new process to fill incomplet genes using MAFFT alignment.
- A more complete command line interface with a sub-command system.
- A configuration file for defined PATH to external tools.
- An easy installation with pypi repository.


## Documentation
The details of installation, workflow and running parameters could be found on the [**documentation**](https://pymlst.readthedocs.io/en/latest/).


## Publications
If you use pyMLST, please cite the following paper:

Bignenet A. et al., Introduction and benchmarking of pyMLST:
open-source software for assessing bacterial clonality using core
genome MLST. 2023 Microbials Genomics, 9(11), 1126.
doi: [10.1099/mgen.0.001126](https://doi.org/10.1099/mgen.0.001126)


PyMLST v1 have been already use to analyse most of clinical bacteria:

 - [*Escherichia coli* and *Klebsiella pneumoniae*](https://doi.org/10.1016/j.cmi.2021.07.022)
 - [*Acinetobacter baumanii*](https://doi.org/10.1038/s41598-023-49268-x)
 - [*Pseudomonas aeruginosa*](https://doi.org/10.1016/j.jhin.2020.06.013)
 - [*Proteus mirabilis*](https://doi.org/10.1093/jac/dkz472)
