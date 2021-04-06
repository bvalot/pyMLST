import os
import re
import sys

from Bio import SeqIO
from decorator import contextmanager

from sqlalchemy import create_engine
from sqlalchemy import and_

from sqlalchemy import MetaData, Table, Column, Integer, Text

from sqlalchemy.sql import select
from sqlalchemy.sql import distinct

from pymlst.common import blat
from pymlst.common.binaries import get_binary_path
from pymlst.common.utils import read_genome, create_logger


@contextmanager
def open_cla(file=None, ref='ref'):
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

    def __init__(self, path):
        self.engine = create_engine('sqlite:///' + path)

        metadata = MetaData()

        self.sequences = Table('sequences', metadata,
                               Column('id', Integer, primary_key=True),
                               Column('sequence', Text, unique=True),
                               Column('gene', Text),
                               Column('allele', Integer))

        self.mlst = Table('mlst', metadata,
                          Column('id', Integer, primary_key=True),
                          Column('st', Integer),
                          Column('gene', Text),
                          Column('allele', Integer))

        metadata.create_all(self.engine)

        self.connection = self.engine.connect()
        self.transaction = self.connection.begin()

    def add_sequence(self, sequence, gene, allele):
        self.connection.execute(
            self.sequences.insert(),
            sequence=sequence, gene=gene, allele=allele)

    def add_mlst(self, st, gene, allele):
        self.connection.execute(
            self.mlst.insert(),
            st=st, gene=gene, allele=allele)

    def get_genes_by_allele(self, allele):
        """Returns all the distinct genes in the database and their sequences for a given allele"""
        return self.connection.execute(
            select([distinct(self.mlst.c.gene), self.sequences.c.sequence])
            .select_from(self.mlst.join(
                self.sequences,
                self.mlst.c.gene == self.sequences.c.gene))
            .where(self.sequences.c.allele == allele)
        ).fetchall()

    def get_allele_by_sequence_and_gene(self, sequence, gene):
        return self.connection.execute(
            select([self.sequences.c.allele])
            .where(and_(
                self.sequences.c.sequence == sequence,
                self.sequences.c.gene == gene))
        ).fetchone()

    def get_strains_by_gene_and_allele(self, gene, allele):
        return self.connection.execute(
            select([self.mlst.c.st])
            .where(and_(
                self.mlst.c.gene == gene,
                self.mlst.c.allele == allele))
        ).fetchall()

    def commit(self):
        self.transaction.commit()

    def rollback(self):
        self.transaction.rollback()

    def close(self):
        self.engine.dispose()


class ClassicalMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseCLA(file)
        self.ref = ref
        self.blat_path = '/usr/bin/'
        self.logger = create_logger()

    def create(self, scheme, alleles):
        # Verify sheme list with fasta files
        header = scheme.readline().rstrip("\n").split("\t")
        if len(header) != len(alleles) + 1:
            raise Exception("The number of genes in sheme don't correspond to the number of fasta file\n"
                            + " ".join(header) + "\n")
        fastas = {}
        for f in alleles:
            name = f.name.split("/")[-1]
            name = name[:name.rfind(".")]
            if name not in header:
                raise Exception("Gene " + name + " not found in sheme\n" + " ".join(header))
            fastas[name] = f

        # load sequence allele
        alleles = {}
        for g, f in fastas.items():
            alleles[g] = set()
            for seq in SeqIO.parse(f, 'fasta'):
                try:
                    # if len(seq.id.split("_")) == 2:
                    #     allele = int(seq.id.split("_")[1])
                    # elif len(seq.id.split("-")) == 2:
                    #     allele = int(seq.id.split("-")[1])
                    # elif g in seq.id:
                    #     allele = int(seq.id.replace(g, ""))
                    # else:
                    #     allele = int(seq.id)
                    match = re.search('[0-9]+$', seq.id)
                    allele = int(match.group(0))

                except Exception:
                    raise Exception("Unable to obtain allele number for the sequence: " + seq.id)
                self.database.add_sequence(str(seq.seq).upper(), g, allele)
                alleles.get(g).add(allele)

        # load MLST sheme
        for line in scheme:
            h = line.rstrip("\n").split("\t")
            st = int(h[0])
            for g, a in zip(header[1:], h[1:]):
                if int(a) not in alleles.get(g):
                    self.logger.info(
                        "Unable to find the allele number " + a + " for gene " + g + "; replace by 0")
                    self.database.add_mlst(st, g, 0)
                else:
                    self.database.add_mlst(st, g, int(a))

        self.logger.info('Database initialized')

    def search_st(self, genome, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 to 1")

        path = get_binary_path('blat')
        if path is None:
            raise Exception('Unable to locate the Blat executable\n')

        tmpfile, tmpout = blat.blat_tmp()

        try:
            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            self.logger.info("Search coregene with BLAT")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage, self.logger)
            self.logger.info("Finish run BLAT, found " + str(len(genes)) + " genes")

            # Search sequence MLST
            seqs = read_genome(genome)
            self.logger.info("Search allele gene to database")
            # print(genes)
            allele = {i: [] for i in coregenes}
            st = {i: set() for i in coregenes}
            for coregene in coregenes:
                if coregene not in genes:
                    allele.get(coregene).append("")
                    continue
                for gene in genes.get(coregene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # verify coverage and correct
                    if gene.coverage != 1:
                        gene.searchCorrect()
                        self.logger.info("Gene " + gene.geneId() + " fill: added")

                    # get sequence
                    sequence = str(gene.getSequence(seq)).upper()

                    # verify complet sequence
                    if len(sequence) != (gene.end - gene.start):
                        self.logger.info("Gene " + gene.geneId() + " removed")
                        continue

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
                        strains = self.database.get_strains_by_gene_and_allele(coregene, res[0])
                        for strain in strains:
                            st.get(coregene).add(strain[0])
                    else:
                        allele.get(gene.geneId()).append("new")

            # if only know allele or not found
            # Seach st
            st_val = []
            if sum([len(i) == 1 and i[0] != "new" for i in allele.values()]) == len(allele):
                tmp = None
                for s in st.values():
                    if s:
                        if tmp is None:
                            tmp = s
                        else:
                            tmp = tmp.intersection(s)
                st_val = list(tmp)

            # print result
            coregenes.sort()
            output.write("Sample\tST\t" + "\t".join(coregenes) + "\n")
            output.write(genome.name + "\t" + ";".join(map(str, st_val)))
            for coregene in coregenes:
                output.write("\t" + ";".join(map(str, allele.get(coregene))))
            output.write("\n")
            self.logger.info("FINISH")
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def __create_coregene(self, tmpfile):
        ref = int(1)
        all_rows = self.database.get_genes_by_allele(ref)
        coregenes = []
        for row in all_rows:
            tmpfile.write('>' + row[0] + "\n" + row[1] + "\n")
            coregenes.append(row[0])
        return coregenes

    def close(self):
        self.database.close()

    def commit(self):
        self.database.commit()

    def rollback(self):
        self.database.rollback()