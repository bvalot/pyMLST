from pymlst.wg_commands.db.model import Base, Mlst, Sequence
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError


class Database:

    def __init__(self, path, create=False):
        self.engine = create_engine('sqlite:///' + path)
        self.session_factory = sessionmaker(bind=self.engine)
        self.session = self.session_factory()  # Retrieves a session from a pool maintained by the Engine

        if create:
            Base.metadata.create_all(self.engine)

    def add_mlst(self, souche, gene, seqid):
        self.session.add(Mlst(souche=souche, gene=gene, seqid=seqid))

    def add_sequence(self, sequence):
        entry = Sequence(sequence=sequence)
        try:
            # res = self.session.execute(
            #     Sequence.__table__.insert().values(sequence=sequence)
            # )
            self.session.add(entry)
            self.session.flush()
        except IntegrityError:  # duplicated sequence
            return -1
        return entry.id

    def concatenate_gene(self, gene):
        seq = self.session.query(Sequence) \
                    .filter_by(sequence=str(gene.seq).upper()) \
                    .first()
        existing_gene = self.session.query(Mlst) \
                            .filter_by(seqid=seq.seqid) \
                            .first()
        existing_gene.id += ';' + gene.id

    def remove_sequences(self, sequences):
        ids = []
        for seq in sequences:
            ids.append(self.session.query(Sequence)
                           .filter_by(sequence=seq)
                           .first()
                           .id)
        self.session.query(Sequence) \
                    .filter(Sequence.id.in_(ids)) \
                    .delete()
        self.session.query(Mlst) \
            .filter(Mlst.seqid.in_(ids)) \
            .delete()

    def commit(self):
        self.session.commit()

    def close(self):
        self.session.close()
        self.engine.dispose()

    def rollback(self):
        self.session.rollback()
