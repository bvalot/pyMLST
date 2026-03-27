# Changelog

## 2.3.0 - 2025-03-27
### Added
- gzip compression compatibility for genome fasta file input (fna.gz, fasta.gz)
### Changed
- Removed gene with ambiguous base (N) as valid genes that could result to :
  - over-estimated cgMLST distance between 2 strains
  - create new alleles on classical MST


## 2.2.2 - 2025-06-13
### Fixed
- Mafft patch by @bvalot in #39


## 2.2.1 - 2025-03-27
### Added
- Add new typing method (pyTyper command) 
  - fimH and Clermont typing for Escherichia coli
  - spa typing for Staphylococcus aureus
- Add informations of species and version in database (classical and cg/wgMLST)


## 2.1.6 - 2024-03-20
### Added
- Pasteur repository for claMLST import
### Changed
- Mafft alignment is performed directly instead of using biopython wrapper
### Fixed
- Correct kma search bug on claMLST


## 2.1.5 - 2023-11-09
Release associated to publish the software on Microbial Genomics.
### Fixed
- Correction of performance for distance calculation


## 2.1.4 - 2023-06-30
### Fixed
- Report error on empty database


## 2.1.3 - 2022-06-13
First release with conda installation available
### Fixed
- Correct kma result reading


## 2.1.2 - 2022-03-25
### Added
- Use raw reads (FASTQ) directly with the kma integration (search2 and add2 command)


## 2.0.1 - 2021-09-28
First release of pyMLST V2
### Changed
- An automatic import database mechanism to initiated cgMLST and MLST databases.
- A new process to fill incomplet genes using MAFFT alignment.
- A more complete command line interface with a sub-command system.
- A configuration file for defined PATH to external tools.
- An easy installation with pypi repository.
