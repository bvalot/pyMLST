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
The entry is a draft genome produced by an assembler or the direct raw data, but also other genomes store in sequence database.

## New version
V2.1:

- Use raw reads (FASTQ) directly with the kma integration (search2 and add2 command)

V2.0:

- An automatic import database mechanism to initiated cgMLST and MLST databases.
- A new process to fill incomplet genes using MAFFT alignment.
- A more complete command line interface with a sub-command system.
- A configuration file for defined PATH to external tools.
- An easy installation with pypi repository.


## Installation
### From pypi repository

```
pip install pymlst
```

### From source

```
virtualenv venv
source venv/bin/activate
make install
make build
```

### Dependancy

PyMLST uses 3 external tools to run alignment:

- Mafft (>=7.307)
```
sudo apt install mafft 
```
- Blat (v35).
You need to compile source or obtaine executable at:
[https://genome.ucsc.edu/FAQ/FAQblat.html](https://genome.ucsc.edu/FAQ/FAQblat.html)
- kma (>=1.3)
You need to compile source from:
[https://bitbucket.org/genomicepidemiology/kma/src/master/](https://bitbucket.org/genomicepidemiology/kma/src/master/)


### Configuration

Configure the executables locations and log level :

```
pyMLST configure --help
Usage: pyMLST configure [OPTIONS]

  Configure executables paths and log level.

Options:
  -b, --blat FILE   Blat executable absolute path.
  -k, --kma FILE    Kma executable absolute path.
  -m, --mafft FILE  Mafft executable absolute path.
  -l, --log [DEBUG|INFO|WARNING|ERROR]
                    Level of logging, default=INFO  
  -r, --reset       Reset the configuration.
  --help            Show this message and exit.
```

## cg/wgMLST Analysis

A complete analysis of wgMLST is performed using a succession of python script.

```
wgMLST --help
Usage: wgMLST [OPTIONS] COMMAND [ARGS]...

  Whole/Core genome MLST analysis.

Options:
  -v, --version  Prints PyMLST version.

Commands:
  add            Add a strain GENOME to the wgMLST DATABASE.
  create         Create a wgMLST DATABASE from a template COREGENE.
  distance       Extract an distance matrix from a wgMLST DATABASE.
  gene           Extract an genes list from a wgMLST DATABASE.
  import         Create a wgMLST DATABASE from an online resource.
  mlst           Extract an MLST table from a wgMLST DATABASE.
  msa            Compute Multiple Sequence Alignment from a wgMLST...
  recombinaison  Search potential gene re-combinations from wgMLST...
  remove         Remove STRAINS or GENES from a wgMLST DATABASE.
  sequence       Extract sequences from a wgMLST DATABASE.
  stats          Extract stats from a wgMLST DATABASE.
  strain         Extract an strains list from a wgMLST DATABASE.
  subgraph       Search group of strains at a DISTANCE threshold.

```

### Import or Create a database

First, you need to create a database containing the scheme to use.

- You can simply download existing cgMLST profiles from <https://www.cgmlst.org/> with the **import** command.
```
wgMLST import -h
Usage: wgMLST import [OPTIONS] DATABASE [SPECIES]...

  Create a wgMLST DATABASE from an online resource.

  The research can be filtered by adding a SPECIES name.

Options:
  -f, --force             Override alrealdy existing DATABASE
  --prompt / --no-prompt  Do not prompt if multiple choices are found,
                          fail instead.
```

- Alternatively, you can **create** one with a scheme.
The scheme is a multi-fasta file containing sequences of genes in nucleotide format.
You can obtained scheme for:
	- Core genome analysis in publications.
	- Whole genome analysis by using annoted genes of a publish genome close to your strains.

```
wgMLST create --help
Usage: wgMLST create [OPTIONS] COREGENE DATABASE

  Create a wgMLST DATABASE from a template COREGENE.

Options:
  -f, --force        Override alrealdy existing DATABASE
  -c, --concatenate  Automatically concatenate genes with duplicated sequences
  -r, --remove       Automatically remove genes with duplicated sequences
```

### Add strains

Next, you need to add your strain iteratively to the database. 
A draft genome can be used (we recommend to use [Spades](http://cab.spbu.ru/software/spades/) for assembly).
You can also add reference genome for comparison.

```
wgMLST add --help
Usage: wgMLST add [OPTIONS] GENOME DATABASE

  Add a strain GENOME to the wgMLST DATABASE.

Options:
  -s, --strain TEXT     Name of the strain (default:genome name)
  -i, --identity FLOAT  Minimum identity to search gene (default=0.95)
  -c, --coverage FLOAT  Minimum coverage to search gene (default=0.9)
```

Alternatively, you can also add strain from raw reads direcly with single or paired FASTQS(.gz) files. 

```
wgMLST add2 --help
Usage: wgMLST add2 [OPTIONS] DATABASE [FASTQS]...

  Add a strain from FASTQS(.gz) reads to the wgMLST DATABASE.

Options:
  -s, --strain TEXT     Name of the strain (default:genome name).
  -i, --identity FLOAT  Minimum identity to search gene (default=0.95).
  -c, --coverage FLOAT  Minimum coverage to search gene (default=0.9).
  -r, --reads INTEGER   Minimum reads coverage to search a gene (default=10).
```

### Export results

When the database is complete, you can export results for futher analysis.

- **Distance:**
A matrix of cgMLST distance can be compute from the database and defined as the number of different alleles between each pair of two strains, omitting the missing data.
The genes used to compute this distance can be filtered with the different options (-m, -k, -d and -V). 

```
wgMLST distance --help
Usage: wgMLST distance [OPTIONS] DATABASE

  Extract an distance matrix from a wgMLST DATABASE.
Options:
  -m, --mincover INTEGER  Minimun number of strain found to keep a gene
                          (default:0)

  -k, --keep              Keep only gene with different allele (omit
                          missing).

  -d, --duplicate         Conserve duplicate gene (default remove).
  -V, --inverse           Keep only gene that do not meet the filter of
                          mincover or keep options.

  -o, --output FILENAME   Export distance to (default=stdout).
```

- **MLST**:
The MLST profiles can be also exported. An formatted version compatible with grapetree can be defined.

```
wgMLST mslt --help
Usage: wgMLST mlst [OPTIONS] DATABASE

  Extract an MLST table from a wgMLST DATABASE.
Options:
  ...
  -f, --form [default|grapetree]  Specify format of output
```


### Export sequences

You can access to allele sequences present in the database and specify a list of genes to export with **-l** option (The gene list can be obtained with the **gene** command).

- **Sequence**:
A simple export of the different sequences

```
wgMLST sequence -h
Usage: wgMLST sequence [OPTIONS] DATABASE

  Extract sequences from a wgMLST DATABASE.

Options:
  -o, --output FILENAME  Output result in fasta format (default:stdout).
  -f, --file FILENAME    File containing list of coregenes to extract
                         (default:all coregenes).
  --reference            Return sequence of the reference instead of strains
                         alleles
```

- **MSA**:
A multialign fasta file with genes concatened. 
The file can be use directly for phylogenetic analysis using maximun likelihood or bayesien approaches.

```
wgMLST msa -h
Usage: wgMLST msa [OPTIONS] DATABASE

  Compute Multiple Sequence Alignment from a wgMLST DATABASE.

Options:
  ...
  -r, --realign          Realign genes with same length (Default:No).
```

## classical MLST Analysis

Furthermore, pyMLST is able to search classical MLST and return alleles number and Sequence Type. 

```
claMLST --help
Usage: claMLST [OPTIONS] COMMAND [ARGS]...

  Classical MLST commands.

Options:
  -v, --version  Prints PyMLST version.
  --help         Show this message and exit.

Commands:
  create  Create a classical MLST DATABASE from a SCHEME csv and ALLELES...
  import  Create a claMLST DATABASE from an online resource.
  search  Search ST number for assembly GENOMES using an mlst DATABASE
```

### Initialise a MLST database

- **Import** :
You can import a MLST resource from <https://pubmlst.org/data/>.

```
claMLST import -h
Usage: claMLST import [OPTIONS] DATABASE [SPECIES]...

  Create a claMLST DATABASE from an online resource.

  The research can be filtered by adding a SPECIES name.

Options:
  --prompt / --no-prompt  Do not prompt if multiple choices are found,
                          fail instead.
  -f, --force        	  Override alrealdy existing DATABASE
  -m, --mlst TEXT         Specify the desired MLST scheme name.

```

- **Create** :
Alternatively, you can create a database with the sequence of alleles and mlst profile of your specie of interest at.

To create database, pyMLST needs the gene name present in the mlst profile header to match the name of the fasta file.
As an example, rpoB gene in the header of mlst profile must match rpoB.fas file. 
You also need to remove additionnal column corresponding to clonal complex in the mlst profile file, if existing.

```
claMLST create --help
Usage: claMLST create [OPTIONS] DATABASE SCHEME ALLELES...

  Create a classical MLST DATABASE from a SCHEME csv and ALLELES files.

Options:
  -f, --force        	  Override alrealdy existing DATABASE
```

### Search MLST profile of a strain

Similarly to wgMLST analysis, you need a draft genome to find the mlst profile.
In case a new allele is present, you can obtain the sequence with the **-f** option. Multiple genomes could be search simultanisly.

```
claMLST search --help
Usage: claMLST search [OPTIONS] DATABASE GENOMES

  Search ST number for assembly GENOMES using an mlst DATABASE

Options:
  -i, --identity FLOAT   Minimum identity to search gene (default=0.9)
  -c, --coverage FLOAT   Minimum coverage to search gene (default=0.9)
  -f, --fasta FILENAME   Write fasta file with gene allele
  -o, --output FILENAME  Write ST search result to (default:stdout)
```

Alternatively, you can also search ST from raw reads direcly with single or paired FASTQS(.gz) files.

```
claMLST search2 --help
Usage: claMLST search2 [OPTIONS] DATABASE [FASTQS]...

  Search ST number from FASTQS(.gz) raw reads using an mlst DATABASE.

Options:
  -i, --identity FLOAT   Minimum identity to search gene (default=0.9).
  -c, --coverage FLOAT   Minimum coverage to search gene (default=0.95).
  -r, --reads INTEGER    Minimum reads coverage to search gene (default=10).
  --paired / --single    Defined type of fastqs files.
  -f, --fasta FILENAME   Write fasta file with gene allele.
  -o, --output FILENAME  Write ST search result to (default:stdout).
```

# Publications
PyMLST v1 have been already use to analyse most of clinical bacteria:

 - [*Escherichia coli* and *Klebsiella pneumoniae*](https://doi.org/10.1016/j.cmi.2021.07.022)
 - *Acinetobacter baumanii* (in review)
 - [*Pseudomonas aeruginosa*](https://doi.org/10.1016/j.jhin.2020.06.013)
 - [*Proteus mirabilis*](https://doi.org/10.1093/jac/dkz472)