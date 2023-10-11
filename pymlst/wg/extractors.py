"""Set of methods to extract different types of results from wgMLST"""
import abc
import importlib
import logging
import click

from abc import ABC
import pandas as pd

from pymlst.common import mafft, exceptions, utils
from pymlst.wg.core import Extractor


def read_gene_list(base, gene_file):
    core = base.get_core_genes()
    if gene_file is None:
        return core
    else:
        select = []
        for g in utils.strip_file(gene_file):
            if g in core:
                select.append(g)
            else:
                logging.debug("Gene {} not found in the database".format(g))
        return select


class SequenceExtractor(Extractor):
    """ Extracts coregene sequences into fasta file."""

    def __init__(self, file=None, reference=False):
        """    
        :param file: Path of the file containing the coregens to extract 
        """
        self.list_file = file
        self.reference = reference

    def extract(self, base, output):
        coregene = read_gene_list(base, self.list_file)
        logging.info("Number of gene to analyse : %s", len(coregene))
        for gene in coregene:
            if self.reference:
                seq = base.get_gene_sequence_reference(gene)
                output.write(">" + gene + "|reference" + "\n")
                output.write(seq + "\n")                
            else:
                seqs = base.get_gene_sequences(gene)
                for seq in seqs:
                    output.write(">" + gene + "|" + str(seq[0]) + " "
                                 + ";".join(seq[1]) + "\n")
                    output.write(seq[2] + "\n")

                    
class MsaExtractor(Extractor):
    """ Compute Multiple Sequence Alignment (MSA) and extracts the aligned sequences. """

    def __init__(self, file=None, realign=False):
        """
        :param file: Path of the file containing the coregens to extract
        :param realign: Realign genes with same length 
        """
                
        self.list_file = file
        self.realign = realign

    def extract(self, base, output):
        coregene = read_gene_list(base, self.list_file)
        if len(coregene) == 0:
            raise exceptions.PyMLSTError('No valid genes selected, verify your genes list')
        strains = base.get_all_strains()
        duplicated = base.get_duplicated_genes()

        sequences = {s: [] for s in strains}
        for index, gene in enumerate(coregene):
            if gene in duplicated:
                logging.info("%s/%s | %s     %s", index + 1, len(coregene), gene, "No: Repeat gene")
                continue
            seqs = base.get_gene_sequences(gene)
            size = set()
            for seq in seqs:
                size.add(len(seq[2]))
            if len(size) == 1 and self.realign is False:
                self.add_sequence_strain(seqs, strains, sequences)
                logging.info("%s/%s | %s     %s", index + 1, len(coregene), gene, "Direct")
            else:
                genes = {str(s[0]): s[2] for s in seqs}
                corrseqs = mafft.align(genes)
                for seq in seqs:
                    seq[2] = corrseqs.get(str(seq[0]))
                self.add_sequence_strain(seqs, strains, sequences)
                logging.info("%s/%s | %s     %s", index + 1, len(coregene), gene, "Align")

        # output align result
        for strain in strains:
            output.write('>' + strain + "\n")
            output.write("\n".join(map(str, sequences.get(strain))) + "\n")

    def add_sequence_strain(self, seqs, strains, sequences):
        """Add a sequence to multi-align, take the first gene in case of repetition"""
        size = 0
        if len(seqs) > 0:
            size = len(seqs[0][2])
        for strain in strains:
            seq = [i[2] for i in seqs if strain in i[1]]
            if len(seq) == 0:
                sequences.get(strain).append('-' * size)
            elif len(seq) == 1:
                sequences.get(strain).append(seq[0])
            else:
                raise exceptions.PyMLSTError(
                    'Repeated genes must be excluded in order to export alignment')


            

class TableExtractor(Extractor):
    """ Extraction of cgMLST distance matrix, MLST profiles, Genes and Strains list from a wgMLST database.  """
    def __init__(self,
                 mincover=0,
                 keep=False,
                 duplicate=False,
                 inverse=False):
        self.mincover = mincover
        self.keep = keep
        self.duplicate = duplicate
        self.inverse = inverse

    @abc.abstractmethod
    def extract(self, base, output):
        pass

    def get_valid_shema(self, base):
        # read samples mlst
        strains = base.get_all_strains()
        # Minimun number of strain
        if self.mincover < 0 or self.mincover > len(strains):
            raise exceptions.PyMLSTError(
                'Mincover must be between 0 and number of strains {}'.format(len(strains)))

        # allgene
        allgene = base.get_core_genes()
        # duplicate gene
        dupli = base.get_duplicated_genes()
        # cover without duplication
        count_souches = base.count_souches_per_gene()
        # Count distinct gene
        diff = base.count_sequences_per_gene()

        # filter coregene that is not sufficient mincover or keep only different or return inverse
        valid_shema = []
        # Test different case for validation
        for gene in allgene:
            valid = []
            if self.keep is True:
                if diff.get(gene, 0) > 1:
                    valid.append(True)
                else:
                    valid.append(False)
            else:
                valid.append(True)
            if count_souches.get(gene, 0) >= self.mincover:
                valid.append(True)
            else:
                valid.append(False)
            if not self.duplicate:
                if gene in dupli:
                    valid.append(False)
                else:
                    valid.append(True)
            else:
                valid.append(True)
            if self.inverse is False:
                if sum(valid) == 3:
                    valid_shema.append(gene)
            else:
                if sum(valid) < 3:
                    valid_shema.append(gene)

        # report
        logging.info("Number of coregene used : %s/%s", len(valid_shema), len(allgene))
        return(valid_shema)

class TableExtractorCommand(click.core.Command):
    """ Options supported by :class:`~pymlst.wg.extractors.TableExtractor`. """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params.insert(0, click.core.Option(('--mincover', '-m'),
            type=click.INT,
            help='Minimum number of strains found to retain a gene (default:0)'))
        self.params.insert(1, click.core.Option(('--keep', '-k'),
            is_flag=True,
            help='Keeps only gene with different alleles (omit missing).'))
        self.params.insert(2, click.core.Option(('--duplicate', '-d'),
            is_flag=True,
            help='Keeps duplicate genes (default remove).'))
        self.params.insert(3, click.core.Option(('--inverse', '-V'),
            is_flag=True,
            help='Keeps only gene that do not ' \
                'match the filter of mincover or keep options.'))

class GeneExtractor(TableExtractor):
    """ Extracts a list of genes from a wgMLST database. """
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        
    def extract(self, base, output):
        valid_schema = super().get_valid_shema(base)
        output.write("\n".join(sorted(valid_schema)) + "\n")

class StatsExtractor(Extractor):
    """ Extracts stats, number of strains, coregenes and sequences from a wgMLST database. """
    def extract(self, base, output):
        output.write("Strains\t" + str(len(base.get_all_strains())) + "\n")
        output.write("Coregenes\t" + str(len(base.get_core_genes())) + "\n")
        output.write("Sequences\t" + str(base.count_sequences()) + "\n")

class StrainExtractor(TableExtractor):
    """ Extracts a list of strains from a wgMLST database. """
    def __init__(self, count=False, **kwargs):
        super().__init__(**kwargs)
        self.count = count
        
    def extract(self, base, output):
        if self.count is False:
            output.write("\n".join(base.get_all_strains()) + "\n")
        else:
            tmp = base.count_genes_per_souche(super().get_valid_shema(base))
            for strain in base.get_all_strains():
                output.write(strain + "\t" + str(tmp.get(strain)) + "\n")

class DistanceExtractor(TableExtractor):
    """ Extracts a distance matrix from a wgMLST database. """    
    def extract(self, base, output):
        if self.duplicate:
            logging.warning("Calculate distance between strains " +
                         "using duplicate genes could reported bad result.")
        strains = base.get_all_strains()
        output.write(str(len(strains)) + "\n")
        distance = base.get_strains_distances(super().get_valid_shema(base))
        for strain in strains:
            output.write(strain + "\t")
            dist = [str(distance.get(strain, {}).get(s2, 0)) for s2 in strains]
            output.write("\t".join(dist) + "\n")

class MlstExtractor(TableExtractor):
    """ Extracts an MLST table from a wgMLST database. """
    def __init__(self, form="default", **kwargs):
        super().__init__(**kwargs)
        self.form = form

    def extract(self, base, output):
        valid_shema = super().get_valid_shema(base)
        strains = base.get_all_strains()
        mlst = base.get_mlst(valid_shema)
        table = pd.DataFrame(columns=["#GeneId"] + strains)
        rows = []
        for gene in valid_shema:
            row = {"#GeneId": gene}
            mlstg = mlst.get(gene, {})
            for strain in strains:
                row[strain] = mlstg.get(strain, None)
            rows.append(row)
        table = pd.concat([table, pd.DataFrame.from_dict(rows)], ignore_index=True)
        table = table.set_index('#GeneId')
        
        if self.form == 'grapetree':
            if self.duplicate:
                logging.warnings("Export grapetree table " +
                             "using duplicate genes is not recommended.")
            table = table.fillna(-1)
            table = table.transpose()
        else:
            table = table.fillna("")
        
        table.to_csv(output, sep='\t')
