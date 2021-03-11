import sys
from contextlib import contextmanager

from Bio import SeqIO

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


class WholeGenomeMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseWG(file)
        self.ref = ref

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

    def close(self):
        self.database.close()

    def commit(self):
        self.database.commit()

    def rollback(self):
        self.database.rollback()
