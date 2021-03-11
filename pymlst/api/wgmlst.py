import os
import sys
import tempfile
from contextlib import contextmanager

from Bio import SeqIO

from pymlst.lib import blat, psl
from pymlst.wg_commands.db.database import DatabaseWG


@contextmanager
def open_wg(file=None, ref='ref'):
    mlst = WholeGenomeMLST(file, ref)
    try:
        yield mlst
    except Exception:
        mlst.rollback()
        raise
    else:
        mlst.commit()
    finally:
        mlst.close()


def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs


class WholeGenomeMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseWG(file)
        self.ref = ref
        self.blat_path = '/usr/bin/'  # TODO: change this

    def create(self, coregene, concatenate=True, remove=True):
        genes = set()
        to_remove = set()

        for gene in SeqIO.parse(coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID: " + gene.id)
            else:
                genes.add(gene.id)

            added, seq_id = self.database.add_sequence(str(gene.seq))

            if not added:
                if concatenate:
                    self.database.concatenate_gene(seq_id, gene.id)
                    sys.stderr.write("Concatenate gene " + gene.id + "\n")
                elif remove:
                    to_remove.add(seq_id)
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")
            else:
                self.database.add_mlst(self.ref, gene.id, seq_id)

        if to_remove:
            self.database.remove_sequences(to_remove)
            sys.stderr.write("Remove duplicate sequence: " + str(len(to_remove)) + "\n")

    def add_strain(self, strain, identity, coverage, genome):
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 and 1")
        path = blat.test_blat_exe(self.blat_path)

        name = strain
        if name is None:
            name = genome.name.split('/')[-1]
        if ";" in name:
            raise Exception("Strain name must not contains special ';'\n")

        tmpfile, tmpout = blat.blat_tmp()

        try:
            # verify that the strain is not already in the database
            if self.database.contains_souche(name):
                raise Exception("Strain is already present in database:\n" + name)

            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            sys.stderr.write("Search coregene with BLAT\n")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage)
            sys.stderr.write("Finish run BLAT, found " + str(len(genes)) + " genes\n")

            # add MLST sequence
            seqs = read_genome(genome)
            bad = 0
            for coregene in coregenes:
                if coregene not in genes:
                    continue
                for gene in genes.get(coregene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # Correct coverage
                    if gene.coverage != 1:
                        if gene.searchPartialCDS(seq, coverage) is False:
                            sys.stderr.write("Gene " + gene.geneId() + " partial: removed\n")
                            bad += 1
                            continue
                        else:
                            sys.stderr.write("Gene " + gene.geneId() + " fill: added\n")

                    # Verify CDS
                    if psl.testCDS(gene.getSequence(seq), False) is False:
                        if gene.searchCorrectCDS(seq, coverage) is False:
                            sys.stderr.write("Gene " + gene.geneId() + " not correct: removed\n")
                            bad += 1
                            continue
                        else:
                            sys.stderr.write("Gene " + gene.geneId() + " correct: added\n")

                    # add sequence and MLST
                    sequence = gene.getSequence(seq)

                    # Insert data in database
                    seqid = self.database.add_sequence(str(sequence))[1]
                    self.database.add_mlst(name, gene.geneId(), seqid)

            sys.stderr.write("Add " + str(len(genes) - bad) + " new MLST gene to database\n")
            sys.stderr.write("FINISH\n")

        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def __create_coregene(self, tmpfile):
        ref_genes = self.database.get_gene_by_souche(self.ref)
        coregenes = []
        for row in ref_genes:
            tmpfile.write('>' + row.gene + "\n" + row.sequence + "\n")
            coregenes.append(row[0])
        return coregenes

    def close(self):
        self.database.close()

    def commit(self):
        self.database.commit()

    def rollback(self):
        self.database.rollback()
