"""Core classes and functions to work with Classical MLST data."""

import logging
import os
import re
import sys

from Bio import SeqIO
from decorator import contextmanager

from sqlalchemy import and_

from sqlalchemy.sql import select
from sqlalchemy.sql import distinct

from pymlst.cla import model
from pymlst.common import blat, utils
from pymlst.common.utils import read_genome, create_logger


@contextmanager
def open_cla(file=None, ref=1):
    """A context manager function to wrap the creation a
       :class:`~pymlst.cla.core.ClassicalMLST` object.

    Context managers allow you to instantiate objects using the ``with`` keyword,
    this way you don't have to manage exceptions and the committing/closing processes yourself.

    :param file: The path to the database file to work with.
    :param ref: The name that will be given to the reference strain in the database.

    :yields: A :class:`~pymlst.cla.core.ClassicalMLST` object.
    """
    mlst = ClassicalMLST(file, ref)
    try:
        yield mlst
    except Exception:
        mlst.rollback()
        raise
    else:
        mlst.commit()
    finally:
        mlst.close()


class DatabaseCLA:
    """A core level class to manipulate the genomic database.

        .. warning:: Shouldn't be instantiated directly,
                     see :class:`~pymlst.cla.core.ClassicalMLST` instead.
    """

    def __init__(self, path):
        """
        :param path: The path to the database file to work with.
        """

        self.__engine = utils.get_updated_engine(path, 'cla')
        self.__connection = self.__engine.connect()
        self.__transaction = self.__connection.begin()

    @property
    def connection(self):
        return self.__connection

    def add_sequence(self, sequence, gene, allele):
        """Adds a new sequence associated to a gene and an allele."""
        self.connection.execute(
            model.sequences.insert(),
            sequence=sequence, gene=gene, allele=allele)

    def add_mlst(self, sequence_typing, gene, allele):
        """Adds a new sequence typing, associated to a gene and an allele."""
        self.connection.execute(
            model.mlst.insert(),
            st=sequence_typing, gene=gene, allele=allele)

    def get_genes_by_allele(self, allele):
        """Returns all the distinct genes in the database and their sequences for a given allele."""
        return self.connection.execute(
            select([distinct(model.mlst.c.gene), model.sequences.c.sequence])
                .select_from(model.mlst.join(
                model.sequences,
                model.mlst.c.gene == model.sequences.c.gene))
                .where(model.sequences.c.allele == allele)
        ).fetchall()

    def get_allele_by_sequence_and_gene(self, sequence, gene):
        """Gets an allele by sequence and gene."""
        return self.connection.execute(
            select([model.sequences.c.allele])
                .where(and_(
                model.sequences.c.sequence == sequence,
                model.sequences.c.gene == gene))
        ).fetchone()

    def get_st_by_gene_and_allele(self, gene, allele):
        """Gets a strain by gene and allele."""
        return self.connection.execute(
            select([model.mlst.c.st])
                .where(and_(
                model.mlst.c.gene == gene,
                model.mlst.c.allele == allele))
        ).fetchall()

    def get_sequence_by_gene_and_allele(self, gene, allele):
        """Gets a sequence by gene and allele."""
        return self.connection.execute(
            select([model.sequences.c.sequence])
                .where(and_(
                model.sequences.c.gene == gene,
                model.sequences.c.allele == allele))
        ).fetchone()

    def commit(self):
        """Commits the modifications."""
        self.__transaction.commit()

    def rollback(self):
        """Rollback the modifications."""
        self.__transaction.rollback()

    def close(self):
        """Closes the database engine."""
        self.__engine.dispose()


class ClassicalMLST:
    """Classical MLST python representation.

        Example of usage::

            open_cla('database.db') as db:
                db.create(open('scheme.txt'), [open('gene1.fasta'),
                                               open('gene2.fasta'),
                                               open('gene3.fasta')])
                db.search_st(open('genome.fasta'))
        """

    def __init__(self, file, ref):
        """
        :param file: The path to the database file to work with.
        :param ref: The name that will be given to the reference strain in the database.
        """
        self.database = DatabaseCLA(file)
        self.ref = ref

        create_logger()

    def create(self, scheme, alleles):
        """Creates a classical MLST database from an MLST profile and a list of alleles.

        :param scheme: The MLST profile
        :param alleles: A list of alleles files.

        The MLST profile should be a **CSV** file respecting the following format:

        .. csv-table:: MLST Profile CSV
            :header: "ST", "gene1", "gene2", "gene3", "..."
            :widths: 5, 5, 5, 5, 5

            1, 1, 1, 1, ...
            2, 3, 3, 2, ...
            3, 1, 2, 1, ...
            4, 1, 1, 3, ...
            ..., ..., ..., ..., ...

        """
        # Verify sheme list with fasta files
        header = scheme.readline().rstrip("\n").split("\t")
        if len(header) != len(alleles) + 1:
            raise Exception("The number of genes in sheme don't "
                            "correspond to the number of fasta file\n"
                            + " ".join(header) + "\n")
        fastas = {}
        for file in alleles:
            name = file.name.split("/")[-1]
            name = name[:name.rfind(".")]
            if name not in header:
                raise Exception("Gene " + name + " not found in sheme\n" + " ".join(header))
            fastas[name] = file

        # load sequence allele
        alleles = {}
        for gene, file in fastas.items():
            alleles[gene] = set()
            for seq in SeqIO.parse(file, 'fasta'):

                try:
                    match = re.search('[0-9]+$', seq.id)
                    allele = int(match.group(0))
                except Exception as err:
                    raise Exception("Unable to obtain allele number for the sequence: " + seq.id) from err

                self.database.add_sequence(str(seq.seq).upper(), gene, allele)
                alleles.get(gene).add(allele)

        # load MLST sheme
        for line in scheme:
            line_content = line.rstrip("\n").split("\t")
            sequence_typing = int(line_content[0])
            for gene, allele in zip(header[1:], line_content[1:]):
                if int(allele) not in alleles.get(gene):
                    logging.info(
                        "Unable to find the allele number %s"
                        " for gene %s; replace by 0", str(allele), gene)
                    self.database.add_mlst(sequence_typing, gene, 0)
                else:
                    self.database.add_mlst(sequence_typing, gene, int(allele))

        logging.info('Database initialized')

    def search_st(self, genome, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        """Search the **Sequence Type** number of a strain.

        :param genome: The strain genome we want to add as a `fasta`_ file.
        :param identity: Sets the minimum identity used by `BLAT`_
                         for sequences research (in percent).
        :param coverage: Sets the minimum accepted coverage for found sequences.
        :param fasta: A file where to export genes alleles results in a fasta format.
        :param output: An output for the sequence type research results.
        """
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 to 1")

        tmpfile, tmpout = blat.blat_tmp()

        try:
            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            logging.info("Search coregene with BLAT")
            genes = blat.run_blat(genome, tmpfile, tmpout, identity, coverage)
            logging.info("Finish run BLAT, found %s genes", str(len(genes)))

            # Search sequence MLST
            seqs = read_genome(genome)
            logging.info("Search allele gene to database")
            # print(genes)
            allele = {i: [] for i in coregenes}
            sequence_type = {i: set() for i in coregenes}
            for coregene in coregenes:
                if coregene not in genes:
                    allele.get(coregene).append("")
                    continue
                for gene in genes.get(coregene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # verify coverage and correct
                    # if gene.coverage != 1:
                    #     gene.searchCorrect()
                    #     logging.info("Gene %s fill: added", gene.gene_id())

                    if gene.coverage == 1:
                        sequence = gene.get_sequence(seq)
                    else:
                        coregene_seq = self.database.get_sequence_by_gene_and_allele(
                            coregene, self.ref)[0]
                        sequence = gene.get_aligned_sequence(seq, coregene_seq)

                    sequence = str(sequence)

                    # get sequence
                    # sequence = str(gene.get_sequence(seq)).upper()

                    # verify complet sequence
                    # if not (sequence and len(sequence) == (gene.end - gene.start)):
                    #     print('Len Seq: {} and gene : {} \n -> {}'.format(len(sequence), gene.end - gene.start, sequence))
                    #     logging.info("Gene %s removed", gene.gene_id())
                    #     #continue
                    # print('OK: {}'.format(sequence))

                    # write fasta file with coregene
                    if fasta is not None:
                        fasta.write(">" + coregene + "\n")
                        fasta.write(sequence + "\n")

                    # search allele
                    res = self.database.get_allele_by_sequence_and_gene(sequence, coregene)
                    if res is not None:
                        allele.get(coregene).append(str(res[0]))
                        # cursor.execute('''SELECT st FROM mlst WHERE gene=? and allele=?''',
                        #                (coregene, row[0]))
                        seq_types = self.database.get_st_by_gene_and_allele(coregene, res[0])
                        for seq_type in seq_types:
                            sequence_type.get(coregene).add(seq_type[0])
                    else:
                        allele.get(gene.gene_id()).append("new")

            # if only know allele or not found
            # Seach st
            st_val = []
            if sum([len(i) == 1 and i[0] != "new" for i in allele.values()]) == len(allele):
                tmp = None
                for st_value in sequence_type.values():
                    if st_value:
                        if tmp is None:
                            tmp = st_value
                        else:
                            tmp = tmp.intersection(st_value)
                st_val = list(tmp)

            # print result
            coregenes.sort()
            output.write("Sample\tST\t" + "\t".join(coregenes) + "\n")
            output.write(genome.name + "\t" + ";".join(map(str, st_val)))
            for coregene in coregenes:
                output.write("\t" + ";".join(map(str, allele.get(coregene))))
            output.write("\n")
            logging.info("FINISH")
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def __create_coregene(self, tmpfile):
        all_rows = self.database.get_genes_by_allele(self.ref)
        coregenes = []
        for row in all_rows:
            tmpfile.write('>' + row[0] + "\n" + row[1] + "\n")
            coregenes.append(row[0])
        return coregenes

    def commit(self):
        """Commits the modifications."""
        self.database.commit()

    def rollback(self):
        """Rollback the modifications."""
        self.database.rollback()

    def close(self):
        """Closes the database engine."""
        self.database.close()
