# pyMLST
python Mlst Local Search Tool

## Purpose
Typing of bacteria is an important task of public health in hospital. 
The use of next generation sequencing for analyse or survey epidimic strain is in expansion.
For this purpose, core or whole genome Multilocus Sequence Typing (cgMLST / wgMLST) has become the new standard.
It was an extension of the traditionnal MLST method that use an short list of hoose keeping genes. 
Here, a large set of gene corresponding to the core or the whole genome is used.
Similarly to MLST, each unique sequence corresponds to an specific allele and 
the combination of allele determines the sequence type (ST) of the strain.

pyMLST have been developped to performed this task. 
In comparaison to other tools, it used an local sqlite database to store allele sequences and mlst profiles. 
This permits to iteratively enlarge the collection of strains to compare. 
The entry is a draft genome produced by assembler, but also other genome store in sequence database.

## Installation

pyMLST need an additionnal python library to run:
- Biopython (>=1.68)
```
sudo apt install python-biopython 
```

pyMLST used 2 external tools to run alignment :
- Mafft (>=7.307)
```
sudo apt install mafft 
```
- Blat (v35).
You need to compiled source or obtained executable at:
[https://genome.ucsc.edu/FAQ/FAQblat.html](https://genome.ucsc.edu/FAQ/FAQblat.html)

## cg/wgMLST Analysis

An complete analysis of wgMLST were performed by using a succession of python script.  

### Create database

First, you need to create a database containing the schema to used. 
The schema is a multi-fasta files containing sequence of genes in nucleotide format.
You can obtained schema for :
- Core genome analysis as described <https://www.cgmlst.org/> or in publications.
- Whole genome analysis by using genes annoted in an publish genome closed to your strains.

```
mlst_create_database.py -h
usage: mlst_create_database.py [options] coregene database

Create a wgMLST database from a template

positional arguments:
  coregene    Coregene fasta file as template of MLST
  database    Sqlite database to stock MLST

optional arguments:
  -h, --help  show this help message and exit
```

### Add strain

Next, you need to add your strain iteratively to the database. 
A draft genome is sufficient (we recommended to used [Spades](http://cab.spbu.ru/software/spades/) for assembly).
You can also add reference genome for comparison.

```
./mlst_add_strain.py -h
usage: mlst_add_strain.py [options] genome database

Add a strain to the wgMLST database

positional arguments:
  genome                Genome of the strain
  database              Sqlite database to stock MLST

optional arguments:
  -h, --help            show this help message and exit
  -s [STRAIN], --strain [STRAIN]
                        Name of the strain (default:genome name)
  -i [IDENTITY], --identity [IDENTITY]
                        Minimun identity to search gene (default=0.95)
  -p [PATH], --path [PATH]
                        Path to BLAT executable (default=/usr/bin)
```

### Export results (profiles / distances)

This script permits to obtained informations on the database and extract result depending of the **-e** option selected.
- **stat**: List the contains in database
- **strain**: List the strains present in the database, 
in combination with **-c** options, you can obtained the number of gene found for each strains.
- **gene**: List the genes present in the database. 
In combination with **-m** option, you can restrict to the gene present in the majority of strains.
- **mlst**: Table containing the mlst profile of each strains. 
To simplify the result, you can export only genes with different alleles with the **-k** option.
- **distance**: Matrix of distances between strains. 
Each distance between a paired of strains are calculated as the number of genes with a different alleles omitting the missing data.
For the calculation, you can used only genes present in a sufficient number of strains with the **-m** option. 

```
./mlst_extract_table.py -h
usage: mlst_extract_table.py [options] database

Extract MLST table from an wgMLST database

positional arguments:
  database              Sqlite database to stock MLST

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT], --output [OUTPUT]
                        Export MLST table to (default=stdout)
  -e {mlst,distance,strain,gene,stat}, --export {mlst,distance,strain,gene,stat}
                        Defined the export format; mlst: MLST table (Default);
                        distance: the distance matrix; strain: the strain
                        list; gene: the gene list; stat: statistics of the
                        database
  -c, --count           In strain mode, count the number of gene present in
                        the database
  -m [MINCOVER], --mincover [MINCOVER]
                        Minimun number of strain found to conserved a gene
                        (default:0)
  -k, --keep            Keep only gene with different allele (omit missing)
  -V, --inverse         Conserved only gene that not pass filter of mincover
                        or keep options
```

### Export sequences

This script gives access to allele sequences present in the database. 
You can specify a list of genes to export with **-l** option. 

You can also report an multialign fasta files with genes concatened using **-a** option. 
The file can be used directly for phylogenetic analysis using maximun likelihood or bayesien approaches.

```
./mlst_extract_sequence.py -h                                     
usage: mlst_extract_sequence.py [options] database

Get sequences from wgMLST database

positional arguments:
  database              Sqlite database to stock MLST

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT], --output [OUTPUT]
                        Output result on fasta format in (default:stdout)
  -l [LISTE], --liste [LISTE]
                        List of coregene to extract (default:all)
  -a, --align           Report a concatened multi-fasta file instead of only
                        gene files (default:No)
  -r, --realign         Realign gene with same length (Default:No)
  -m [MINCOVER], --mincover [MINCOVER]
                        Minimun number of strain found to conserved a coregene
                        (default:1)
  -p [PATH], --path [PATH]
                        Path to mafft executable (default=/usr/bin)
```

## classical MLST Analysis

Additionnaly, pyMLST is able to search classical MLST and return alleles number and Sequence Type. 

### Creation of MLST database

You need to download list of alleles and mlst profile of your specie of interest at <https://pubmlst.org/data/>.

To create database, pyMLST need that the gene name present in the mlst profile header correspond to the name of the fasta file.
As example, rpoB gene indicated in the header of mlst profile must correspond to rpoB.fas file. 
You also need to remove additionnal column corresponding to clonal complex in the mlst profile file, if present.

```
./claMLST_create_database.py -h
usage: claMLST_create_database.py [options] database shema alleles

Create a classical MLST database from a shema

positional arguments:
  database    Sqlite database to stock MLST
  shema       Tabular file containing MLST shema
  alleles     Fasta files containing alleles for each MLST genes, One file per
              genes

optional arguments:
  -h, --help  show this help message and exit
```

### Search MLST profile of a strain

Similarly to wgMLST analysis, you need a draft genome to find the mlst profile.
In case of a new allele is present, you can obtained the sequence with the **-f** option.

```
./claMLST_search_ST.py -h
usage: claMLST_search_ST.py [options] genome database

Search ST number for a strain

positional arguments:
  genome                Genome of the strain
  database              Sqlite database containing MLST shema

optional arguments:
  -h, --help            show this help message and exit
  -i [IDENTITY], --identity [IDENTITY]
                        Minimun identity to search gene (default=0.9)
  -f FASTA, --fasta FASTA
                        Write fasta file with gene allele
  -p [PATH], --path [PATH]
                        Path to BLAT executable (default=/usr/bin)
  -o OUTPUT, --output OUTPUT
                        Write ST search result to (default=stdout)
```
