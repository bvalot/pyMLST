# PyMLST
A Python Mlst Local Search Tool.

## Purpose
Typing bacteria is an important public health task in hospital. 
The use of next generation sequencing to analyse or survey epidimic strain is rising.
For this purpose, core or whole genome Multilocus Sequence Typing (cgMLST / wgMLST) has become the new standard.
It is an extension of the traditionnal MLST method that uses a short list of housekeeping genes. 
Here, a large set of gene corresponding to the core or the whole genome is used.
Similarly to MLST, each unique sequence matches a specific allele and 
the combination of allele determines the sequence type (ST) of the strain.

pyMLST have been developped to performed this task. 
In comparaison to other tools, it uses a local sqlite database to store allele sequences and mlst profiles. 
This allows to iteratively enlarge the collection of strains to compare. 
The entry is a draft genome produced by an assembler, but also other genomes store in sequence database.

## Installation

```
pip install pymlst
```

PyMLST uses 2 external tools to run alignment:
- Mafft (>=7.307)
```
sudo apt install mafft 
```
- Blat (v35).
You need to compile source or obtaine executable at:
[https://genome.ucsc.edu/FAQ/FAQblat.html](https://genome.ucsc.edu/FAQ/FAQblat.html)
  
Configure the executables locations :

```
pyMLST configure --help
Usage: pyMLST configure [OPTIONS]

  Configure executables paths.

Options:
  -b, --blat FILE   Blat executable absolute path.
  -m, --mafft FILE  Mafft executable absolute path.
  -r, --reset       Reset the configuration.
  --help            Show this message and exit.
```

## cg/wgMLST Analysis

A complete analysis of wgMLST is performed using a succession of python script.

```
wgMLST --help
Usage: wgMLST [OPTIONS] COMMAND [ARGS]...

  Whole/Core genome MLST commands.

Options:
  -v, --version  Prints PyMLST version.
  -d, --debug    Sets the debug mode ON.
  --help         Show this message and exit.

Commands:
  add_strain          Add a strain GENOME to the wgMLST DATABASE.
  create_db           Create a wgMLST DATABASE from a template COREGENE.
  extract_sequence    Extract sequences from a wgMLST DATABASE.
  extract_table       Extract an MLST table from a wgMLST DATABASE.
  find_recombinaison  Search potential gene re-combinations from wgMLST...
  find_subgraph       Search group of strains at a DISTANCE threshold.
  import              Create a wgMLST DATABASE from an online resource.
  remove_gene         Remove GENES from a wgMLST DATABASE.
  remove_strain       Remove STRAINS from a wgMLST DATABASE.

```

### Create database

First, you need to create a database containing the scheme to use. 
The scheme is a multi-fasta file containing sequences of genes in nucleotide format.
You can obtained scheme for:
- Core genome analysis as described <https://www.cgmlst.org/> or in publications.
- Whole genome analysis by using annoted genes of a publish genome close to your strains.

```
wgMLST create_db --help
Usage: wgMLST create_db [OPTIONS] COREGENE DATABASE

  Create a wgMLST DATABASE from a template COREGENE.

Options:
  -c, --concatenate  Automatically concatenate genes with duplicated sequences
  -r, --remove       Automatically remove genes with duplicated sequences
  --help             Show this message and exit.
```

### Add strain

Next, you need to add your strain iteratively to the database. 
A draft genome can be used (we recommende to use [Spades](http://cab.spbu.ru/software/spades/) for assembly).
You can also add reference genome for comparison.

```
wgMLST add_strain --help
Usage: wgMLST add_strain [OPTIONS] GENOME DATABASE

  Add a strain GENOME to the wgMLST DATABASE.

Options:
  -s, --strain TEXT     Name of the strain (default:genome name)
  -i, --identity FLOAT  Minimum identity to search gene (default=0.95)
  -c, --coverage FLOAT  Minimum coverage to search gene (default=0.9)
  --help                Show this message and exit.
```

### Export results (profiles / distances)

This script allows to obtain informations on the database and extract result depending on the **-e** option selected.
- **stat**: List the content in database.
- **strain**: List the strains present in the database, 
in combination with **-c** option, you can obtained the number of gene found for each strain.
- **gene**: List the genes present in the database. 
In combination with **-m** option, you can restrict to genes present in the majority of strains.
- **mlst**: Table containing the mlst profile of each strain. 
To simplify the result, you can export only genes with different alleles with the **-k** option.
- **distance**: Matrix of distances between strains. 
Each distance between a paire of strains is calculated as the number of genes with a different alleles, omitting the missing data.
For the calculation, you can use only genes present in a sufficient number of strains with the **-m** option. 

```
wgMLST extract_table --help
Usage: wgMLST extract_table [OPTIONS] DATABASE

  Extract an MLST table from a wgMLST DATABASE.

Options:
  -o, --output FILENAME           Export MLST table to (default=stdout)
  -e, --export [strain|gene|distance|mlst|grapetree|stat]
                                  Defined the export format
  -c, --count                     In strain mode, count the number of gene
                                  present in the database

  -m, --mincover INTEGER          Minimun number of strain found to keep a
                                  gene (default:0)

  -k, --keep                      Keep only gene with different allele (omit
                                  missing)

  -d, --duplicate                 Conserve duplicate gene (default remove)
  -V, --inverse                   Keep only gene that do not meet the filter
                                  of mincover or keep options

  --help                          Show this message and exit.
```

### Export sequences

This script gives access to allele sequences present in the database. 
You can specify a list of genes to export with **-l** option. 

You can also report a multialign fasta file with genes concatened using **-a** option. 
The file can be use directly for phylogenetic analysis using maximun likelihood or bayesien approaches.

```
wgMLST extract_sequence --help
Usage: wgMLST extract_sequence [OPTIONS] DATABASE

  Extract sequences from a wgMLST DATABASE.

Options:
  -o, --output FILENAME   Output result in fasta format (default:stdout)
  -l, --list FILENAME     List of coregenes to extract (default:all)
  -a, --align             Report a concatened multi-fasta file instead of only
                          gene files (default:No)

  -r, --realign           Realign genes with same length (Default:No)
  -m, --mincover INTEGER  Minimun number of strain found to keep a coregene
                          (default:1)

  --help                  Show this message and exit.
```

## classical MLST Analysis

Furthermore, pyMLST is able to search classical MLST and return alleles number and Sequence Type. 

```
claMLST --help
Usage: claMLST [OPTIONS] COMMAND [ARGS]...

  Classical MLST commands.

Options:
  -v, --version  Prints PyMLST version.
  -d, --debug    Sets the debug mode ON.
  --help         Show this message and exit.

Commands:
  create_db  Create a classical MLST DATABASE from a SCHEME csv and ALLELES...
  import     Create a claMLST DATABASE from an online resource.
  search_ST  Search ST number for an assembly GENOME using an mlst DATABASE
```

### Creation of MLST database

You need to download list of alleles and mlst profile of your specie of interest at <https://pubmlst.org/data/>.

To create database, pyMLST needs the gene name present in the mlst profile header to match the name of the fasta file.
As an example, rpoB gene in the header of mlst profile must match rpoB.fas file. 
You also need to remove additionnal column corresponding to clonal complex in the mlst profile file, if existing.

```
claMLST create_db --help
Usage: claMLST create_db [OPTIONS] DATABASE SCHEME ALLELES...

  Create a classical MLST DATABASE from a SCHEME csv and ALLELES files.

Options:
  --help  Show this message and exit.
```

### Search MLST profile of a strain

Similarly to wgMLST analysis, you need a draft genome to find the mlst profile.
In case a new allele is present, you can obtain the sequence with the **-f** option.

```
claMLST search_ST --help
Usage: claMLST search_ST [OPTIONS] GENOME DATABASE

  Search ST number for an assembly GENOME using an mlst DATABASE

Options:
  -i, --identity FLOAT   Minimum identity to search gene (default=0.9)
  -c, --coverage FLOAT   Minimum coverage to search gene (default=0.9)
  -f, --fasta FILENAME   Write fasta file with gene allele
  -o, --output FILENAME  Write ST search result to (default:stdout)
  --help                 Show this message and exit.
```

## Importing online resources

In addition of using local data, you can use online resources to initialize **claMLST**
and **wgMLST** databases.

Simply use :
> wgMLST import database.db

Or : 
> claMLST import database.db

And then pick the species you wish to download from the dropdown menu. The used resources are
<https://www.cgmlst.org/> for MLST data and <https://www.pubmlst.org/> for wg/cgMLST.