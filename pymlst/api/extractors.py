import abc
import importlib

import os

import subprocess

import tempfile
from abc import ABC

from pymlst.api.core import Extractor


def run_mafft(path, tmpfile):
    command = [path, '--quiet', tmpfile.name]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    genes = {}
    ids = None
    seq = ""
    for line in iter(proc.stdout.readline, ''):
        if ids is None and line[0] == '>':
            ids = line.lstrip('>').rstrip("\n")
        elif ids is not None and line[0] == '>':
            genes[int(ids)] = seq.upper()
            ids = line.lstrip('>').rstrip("\n")
            seq = ""
        elif ids is not None:
            seq += line.rstrip("\n")
        else:
            raise Exception("A problem occurred while running mafft" + str(line))
    if seq != "":
        genes[int(ids)] = seq.upper()

    return genes


def write_tmp_seqs(tmpfile, seqs):
    tmp = open(tmpfile.name, 'w+t')
    for s in seqs:
        tmp.write(">"+str(s[0])+"\n"+s[2]+"\n")
    tmp.close()


def add_sequence_strain(seqs, strains, sequences):
    """Add a sequence to multi-align, take the first gene in case of repetition"""
    size = 0
    if len(seqs) > 0:
        size = len(seqs[0][2])
    for s in strains:
        seq = [i[2] for i in seqs if s in i[1]]
        if len(seq) == 0:
            sequences.get(s).append('-' * size)
        elif len(seq) == 1:
            sequences.get(s).append(seq[0])
        else:
            raise Exception("repeated gene must be excluded in order to align export\n")


class SequenceExtractor(Extractor):

    def __init__(self, list=None, align=False, realign=False, mincover=1):
        self.list = list
        self.align = align
        self.realign = realign
        self.mincover = mincover
        self.mafft_path = '/usr/bin/mafft'

    def extract(self, base, ref, output, logger):
        tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
        tmpfile.close()
        try:
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

            logger.info("Number of gene to analyse : " + str(len(coregene)) + "\n")

            if self.align is False:
                # no multialign
                for g in coregene:
                    seqs = base.get_gene_sequences(g, ref)
                    for seq in seqs:
                        output.write(">" + g + "|" + str(seq[0]) + " "
                                     + ";".join(seq[1]) + "\n")
                        output.write(seq[2] + "\n")
            else:
                # multialign
                # search duplicate
                dupli = base.get_duplicated_genes(ref)
                for d in dupli:
                    logger.info('Duplicated: ', d)

                sequences = {s: [] for s in strains}
                for i, g in enumerate(coregene):
                    logger.info(str(i + 1) + "/" + str(len(coregene)) + " | " + g + "     ")
                    if g in dupli:
                        logger.info("No: Repeat gene\n")
                        continue
                    seqs = base.get_gene_sequences(g, ref)
                    size = set()
                    for seq in seqs:
                        size.add(len(seq[2]))
                    if len(size) == 1 and self.realign is False:
                        logger.info("Direct")
                        add_sequence_strain(seqs, strains, sequences)
                    else:
                        logger.info("Align")
                        write_tmp_seqs(tmpfile, seqs)
                        corrseqs = run_mafft(self.mafft_path, tmpfile)
                        for seq in seqs:
                            seq[2] = corrseqs.get(seq[0])
                        add_sequence_strain(seqs, strains, sequences)
                    logger.info("\n")

                # output align result
                for s in strains:
                    output.write('>' + s + "\n")
                    output.write("\n".join(sequences.get(s)) + "\n")
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)


class TableExtractor(Extractor):

    def __init__(self, export='mlst', count=False, mincover=0, keep=False, duplicate=True, inverse=False):
        self.export = export
        self.count = count
        self.mincover = mincover
        self.keep = keep
        self.duplicate = duplicate
        self.inverse = inverse

    def extract(self, base, ref, output, logger):
        # read samples mlst
        strains = base.get_all_strains(ref)

        # Minimun number of strain
        if self.mincover < 0 or self.mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : " + str(len(strains)))

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
        for g in allgene:
            valid = []
            if self.keep is True:
                if diff.get(g, 0) > 1:
                    valid.append(True)
                else:
                    valid.append(False)
            else:
                valid.append(True)
            if count_souches.get(g, 0) >= self.mincover:
                valid.append(True)
            else:
                valid.append(False)
            if self.duplicate:
                if g in dupli:
                    valid.append(False)
                else:
                    valid.append(True)
            else:
                valid.append(True)
            if self.inverse is False:
                if sum(valid) == 3:
                    valid_shema.append(g)
            else:
                if sum(valid) < 3:
                    valid_shema.append(g)

        # report
        logger.info("Number of coregene used : " + str(len(valid_shema)) +
                    "/" + str(len(allgene)) + "\n")

        exporter = ExportType.get_type(self.export)

        if exporter is None:
            raise Exception('Unknown export type: ', self.export)

        data = ExportData(valid_shema, strains, allgene, ref, self.count, self.duplicate)
        exporter.export(data, base, output, logger)


class ExportType(ABC):
    @classmethod
    def list_types(cls):
        importlib.import_module('pymlst.api.export_types')
        return [exp.name() for exp in cls.__subclasses__()]

    @classmethod
    def get_type(cls, type_name):
        importlib.import_module('pymlst.api.export_types')
        for type_cls in cls.__subclasses__():
            if type_cls.name() is type_name:
                return type_cls()
        return None

    @abc.abstractmethod
    def export(self, data, base, output, logger):
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
