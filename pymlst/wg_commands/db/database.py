from pymlst.wg_commands.db.model import Base, Mlst, Sequence
from sqlalchemy import create_engine
from sqlalchemy import and_
from sqlalchemy import func
from sqlalchemy import MetaData, Table, Column, Integer, Text, VARBINARY
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select, exists
from sqlalchemy.sql import distinct
from sqlalchemy.sql.operators import in_op as in_
from pymlst.lib.benchmark import benchmark


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
                          Column('souche', Text, index=True),
                          Column('gene', Text, index=True),
                          Column('seqid', Integer, index=True))

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

    def remove_sequences(self, ids):
        """Removes sequences and their associated genes"""
        self.connection.execute(
            self.sequences.delete()
                .where(in_(self.sequences.c.id, ids))
        )
        self.connection.execute(
            self.mlst.delete()
                .where(in_(self.mlst.c.seqid, ids))
        )

    def get_gene_by_souche(self, souche):
        return self.connection.execute(
            select([self.mlst.c.gene, self.sequences.c.sequence])
                .where(and_(
                self.mlst.c.souche == souche,
                self.mlst.c.seqid == self.sequences.c.id
            ))
        ).fetchall()

    def contains_souche(self, souche):
        return self.connection.execute(
            select([self.mlst.c.id])
                .where(self.mlst.c.souche == souche)
                .limit(1)
        ).fetchone() is not None

    @benchmark
    def get_gene_sequences(self, gene, souche):
        res = self.connection.execute(
            select([self.mlst.c.seqid,
                    func.group_concat(self.mlst.c.souche, ';'),
                    self.sequences.c.sequence])
                .where(and_(
                self.mlst.c.seqid == self.sequences.c.id,
                self.mlst.c.souche != souche,
                self.mlst.c.gene == gene))
                .group_by(self.mlst.c.seqid)
        ).fetchall()
        seqs = []
        for seq in res:
            tmp = seq[1].split(";")
            tmp.sort()
            seqs.append([seq[0], tmp, seq[2]])
        return seqs

    def get_different_souches(self, souche):
        return self.connection.execute(
            select([self.mlst.c.souche])
            .where(self.mlst.c.souche != souche)
            .distinct()
        ).fetchall()

    def get_genes_coverages(self, ref):
        return self.connection.execute(
            select([self.mlst.c.gene,
                    func.count(distinct(self.mlst.c.souche))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()

    @benchmark
    def get_duplicated_genes(self, ref):
        m_alias = self.mlst.alias()

        exist_sub = select([self.mlst]) \
            .where(and_(
                self.mlst.c.souche == m_alias.c.souche,
                self.mlst.c.gene == m_alias.c.gene,
                self.mlst.c.id != m_alias.c.id))

        res = self.connection.execute(
            select([self.mlst.c.gene])
            .where(and_(
                exists(exist_sub),
                self.mlst.c.souche != ref
            ))
            .group_by(self.mlst.c.gene)
        ).fetchall()

        return set([row[0] for row in res])

    def commit(self):
        self.transaction.commit()

    def rollback(self):
        self.transaction.rollback()

    def close(self):
        self.engine.dispose()
