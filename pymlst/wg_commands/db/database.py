from pymlst.wg_commands.db.model import Base, Mlst, Sequence
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


class Database:

    def __init__(self, path, create=False):
        self.engine = create_engine('sqlite:///' + path)
        self.session_factory = sessionmaker(bind=self.engine)
        self.session = self.session_factory()  # Retrieves a session from a pool maintained by the Engine

        if create:
            Base.metadata.create_all(self.engine)

    def add_mlst(self, souche, gene, seqid):
        """Adds an MLST gene bound to an existing sequence"""
        self.session.add(Mlst(souche=souche, gene=gene, seqid=seqid))

    def add_sequence(self, sequence):
        """Adds a sequence if it doesn't already exist"""
        existing = self.session.query(Sequence) \
                       .filter(Sequence.sequence == sequence) \
                       .first()

        if existing is not None:
            return False, existing.id

        entry = Sequence(sequence=sequence)
        self.session.add(entry)
        self.session.flush()

        return True, entry.id

    def concatenate_gene(self, seq_id, gene_name):
        """Associates a new gene to an existing sequence using concatenation"""
        existing_gene = self.session.query(Mlst) \
                            .filter_by(seqid=seq_id) \
                            .first()
        existing_gene.gene += ';' + gene_name

    def remove_sequences(self, ids):
        """Removes sequences and their associated genes"""
        self.session.query(Sequence) \
                    .filter(Sequence.id.in_(ids)) \
                    .delete(synchronize_session=False)
        self.session.query(Mlst) \
                    .filter(Mlst.seqid.in_(ids)) \
                    .delete(synchronize_session=False)

    def commit(self):
        self.session.commit()

    def close(self):
        self.session.close()
        self.engine.dispose()

    def rollback(self):
        self.session.rollback()
