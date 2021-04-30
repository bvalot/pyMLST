import abc
import importlib
import logging

from abc import ABC

from pymlst.common import mafft
from pymlst.common.binaries import get_binary_path
from pymlst.wg.core import Extractor


def add_sequence_strain(seqs, strains, sequences):
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
            raise Exception("repeated gene must be excluded in order to align export\n")


class SequenceExtractor(Extractor):

    def __init__(self, gene_list=None, align=False, realign=False, mincover=1):
        self.list = gene_list
        self.align = align
        self.realign = realign
        self.mincover = mincover

    def extract(self, base, ref, output):
        # Minimun number of strain
        strains = base.get_all_strains(ref)
        if self.mincover < 1 or self.mincover > len(strains):
            raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))

        #  Coregene
        coregene = []
        if self.list is not None:
            coregene = [l.rstrip("\n") for l in iter(self.list.readline, '')]
        else:
            coregene = [l[0] for l in base.get_genes_coverages(ref) if l[1] >= self.mincover]

        logging.info("Number of gene to analyse : %s", len(coregene))

        if self.align is False:
            # no multialign
            for gene in coregene:
                seqs = base.get_gene_sequences(gene, ref)
                for seq in seqs:
                    output.write(">" + gene + "|" + str(seq[0]) + " "
                                 + ";".join(seq[1]) + "\n")
                    output.write(seq[2] + "\n")
        else:
            # multialign
            # search duplicate
            duplicated = base.get_duplicated_genes(ref)
            for dupli in duplicated:
                logging.info('Duplicated: ', dupli)

            sequences = {s: [] for s in strains}
            for index, gene in enumerate(coregene):
                logging.info("%s/%s | %s     ", index + 1, len(coregene), gene)
                if gene in duplicated:
                    logging.info("No: Repeat gene\n")
                    continue
                seqs = base.get_gene_sequences(gene, ref)
                size = set()
                for seq in seqs:
                    size.add(len(seq[2]))
                if len(size) == 1 and self.realign is False:
                    logging.info("Direct")
                    add_sequence_strain(seqs, strains, sequences)
                else:
                    logging.info("Align")
                    #write_tmp_seqs(tmpfile, seqs)
                    mafft_path = get_binary_path('mafft')
                    if mafft_path is None:
                        raise Exception('Unable to locate the Mafft executable\n')

                    genes = {str(s[0]): s[2] for s in seqs}
                    corrseqs = mafft.align(genes)
                    for seq in seqs:
                        seq[2] = corrseqs.get(seq[0])
                    add_sequence_strain(seqs, strains, sequences)
                logging.info("\n")

            # output align result
            for strain in strains:
                output.write('>' + strain + "\n")
                output.write("\n".join(sequences.get(strain)) + "\n")


class TableExtractor(Extractor):

    def __init__(self, export='mlst', count=False, mincover=0, keep=False, duplicate=False, inverse=False):
        self.export = export
        self.count = count
        self.mincover = mincover
        self.keep = keep
        self.duplicate = duplicate
        self.inverse = inverse

    def extract(self, base, ref, output):
        # read samples mlst
        strains = base.get_all_strains(ref)

        # Minimun number of strain
        if self.mincover < 0 or self.mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : "
                            + str(len(strains)))

        # allgene
        allgene = base.get_all_genes(ref)

        # duplicate gene
        dupli = base.get_duplicated_genes(ref)

        # cover without duplication
        count_souches = base.count_souches_per_gene(ref)

        # Count distinct gene
        diff = base.count_sequences_per_gene(ref)

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

        exporter = ExportType.get_type(self.export)

        if exporter is None:
            raise Exception('Unknown export type: ', self.export)

        data = ExportData(valid_shema, strains, allgene, ref, self.count, self.duplicate)
        exporter.export(data, base, output)


class ExportType(ABC):
    @classmethod
    def list_types(cls):
        importlib.import_module('pymlst.wg.export_types')
        return [exp.name() for exp in cls.__subclasses__()]

    @classmethod
    def get_type(cls, type_name):
        importlib.import_module('pymlst.wg.export_types')
        for type_cls in cls.__subclasses__():
            if type_cls.name() is type_name:
                return type_cls()
        return None

    @abc.abstractmethod
    def export(self, data, base, output):
        pass

    @staticmethod
    @abc.abstractmethod
    def name():
        return 'undefined'


class ExportData:
    def __init__(self, valid_schema, strains, all_genes, ref, count, duplicate):
        self.valid_schema = valid_schema
        self.strains = strains
        self.all_genes = all_genes
        self.ref = ref
        self.count = count
        self.duplicate = duplicate
