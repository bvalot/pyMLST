from pymlst.wg_commands.db.model import Base, Mlst, Sequence
from sqlalchemy import create_engine
from sqlalchemy import and_
from sqlalchemy import func
from sqlalchemy import MetaData, Table, Column, Integer, Text, VARBINARY
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select
import hashlib


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
        existing = self.session.query(Sequence.id) \
            .filter(Sequence.sequence == sequence) \
            .first()

        if existing is not None:
            return False, existing.id

        new_seq = Sequence(sequence=sequence)
        self.session.add(new_seq)
        self.session.flush()

        return True, new_seq.id

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

    def get_gene_by_souche(self, souche):
        return self.session.query(Mlst.gene, Sequence.sequence) \
            .filter(and_(Mlst.souche == souche, Mlst.seqid == Sequence.id)) \
            .all()

    def contains_souche(self, souche):
        return self.session.query(Mlst) \
                   .filter(Mlst.souche == souche) \
                   .first() is not None

    def get_gene_sequences(self, gene, souche):
        return self.session.query(
            Mlst.seqid,
            func.group_concat(Mlst.souche, ';'),
            Sequence.sequence
        ).filter(and_(
            Mlst.seqid == Sequence.id,
            Mlst.souche != souche,
            Mlst.gene == gene
        )).group_by(Mlst.seqid).all()

    def commit(self):
        self.session.commit()

    def close(self):
        self.session.close()
        self.engine.dispose()

    def rollback(self):
        self.session.rollback()


class DatabaseCore:

    def __init__(self, path):
        self.engine = create_engine('sqlite:///' + path)

        metadata = MetaData()

        self.sequences = Table('sequences', metadata,
                               Column('id', Integer, primary_key=True),
                               Column('sequence', Text, unique=True))

        self.mlst = Table('mlst', metadata,
                          Column('id', Integer, primary_key=True),
                          Column('souche', Text),
                          Column('gene', Text),
                          Column('seqid', Integer))

        metadata.create_all(self.engine)

        self.connection = self.engine.connect()
        self.transaction = self.connection.begin()

    def add_mlst(self, souche, gene, seqid):
        """Adds an MLST gene bound to an existing sequence"""
        self.connection.execute(
            self.mlst.insert(),
            souche=souche, gene=gene, seqid=seqid)

    def add_sequence(self, sequence):
        """Adds a sequence if it doesn't already exist"""
        existing = self.connection.execute(
            select([self.sequences.c.id])
            .where(self.sequences.c.sequence == sequence)
        ).fetchone()

        if existing is not None:
            return False, existing.id

        res = self.connection.execute(
            self.sequences.insert(),
            sequence=sequence)

        return True, res.inserted_primary_key[0]

    def concatenate_gene(self, seq_id, gene_name):
        """Associates a new gene to an existing sequence using concatenation"""
        self.connection.execute(
            self.mlst.update()
                .values(gene=self.mlst.c.gene + ';' + gene_name)
                .where(self.mlst.c.seqid == seq_id))

    def commit(self):
        self.transaction.commit()

    def rollback(self):
        self.transaction.rollback()

    def close(self):
        self.engine.dispose()
