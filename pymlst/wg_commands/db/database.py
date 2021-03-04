from sqlalchemy.sql.functions import count

from pymlst.wg_commands.db.model import Base, Mlst, Sequence
from sqlalchemy import create_engine, text, bindparam, not_
from sqlalchemy import and_
from sqlalchemy import func
from sqlalchemy import MetaData, Table, Column, Integer, Text, ForeignKey
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
                          Column('seqid', Integer, ForeignKey(self.sequences.c.id)))

        metadata.create_all(self.engine)

        self.connection = self.engine.connect()
        self.transaction = self.connection.begin()

        self.cached_queries = {}

    def get_cached_query(self, name, query_supplier):
        if name in self.cached_queries:
            return self.cached_queries[name]
        query = query_supplier()
        self.cached_queries[name] = query
        return query

    def add_mlst(self, souche, gene, seqid):
        """Adds an MLST gene bound to an existing sequence"""
        self.connection.execute(
            self.mlst.insert(),
            souche=souche, gene=gene, seqid=seqid)

    def add_sequence(self, sequence):
        """Adds a sequence if it doesn't already exist"""
        query = self.get_cached_query(
            'add_sequence',
            lambda:
            select([self.sequences.c.id])
                .where(self.sequences.c.sequence == bindparam('sequence')))

        existing = self.connection.execute(
            query,
            sequence=sequence
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
            .where(in_(self.sequences.c.id, ids)))
        self.connection.execute(
            self.mlst.delete()
            .where(in_(self.mlst.c.seqid, ids)))

    def remove_genes(self, genes):
        deleted = []
        for gene in genes:
            seqids = self.connection.execute(
                select([self.mlst.c.seqid])
                .where(self.mlst.c.gene == gene)
            ).fetchall()
            if len(seqids) == 0:
                continue
            self.connection.execute(
                self.mlst.delete()
                .where(self.mlst.c.gene == gene))
            for seqid in seqids:
                self.connection.execute(
                    self.sequences.delete()
                        .where(and_(
                            not_(exists(
                                select([self.mlst.c.id])
                                .where(self.mlst.c.seqid == self.sequences.c.id))),
                            self.sequences.c.id == seqid[0])))
            deleted.append(gene)
        return deleted

    def get_gene_by_souche(self, souche):
        return self.connection.execute(
            select([self.mlst.c.gene, self.sequences.c.sequence])
                .where(and_(
                self.mlst.c.souche == souche,
                self.mlst.c.seqid == self.sequences.c.id))
        ).fetchall()

    def contains_souche(self, souche):
        return self.connection.execute(
            select([self.mlst.c.id])
                .where(self.mlst.c.souche == souche)
                .limit(1)
        ).fetchone() is not None

    def get_gene_sequences(self, gene, souche):
        query = self.get_cached_query(
            'get_gene_sequences',
            lambda:
            select([self.mlst.c.seqid,
                    func.group_concat(self.mlst.c.souche, bindparam('separator')),
                    self.sequences.c.sequence])
            .select_from(
                self.mlst.join(self.sequences))
            .where(and_(
                self.mlst.c.souche != bindparam('souche'),
                self.mlst.c.gene == bindparam('gene')))
            .group_by(self.mlst.c.seqid))

        res = self.connection.execute(
            query,
            separator=';',
            souche=souche,
            gene=gene
        ).fetchall()

        # res = self.connection.execute(
        #     '''SELECT mlst.seqid, group_concat(mlst.souche, ";") AS group_concat_1, sequences.sequence
        #     FROM mlst JOIN sequences ON sequences.id = mlst.seqid
        #     WHERE mlst.souche != ? AND mlst.gene = ? GROUP BY mlst.seqid''', (souche, gene)
        # )

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

    def get_all_strains(self, ref):
        res = self.connection.execute(
            select([distinct(self.mlst.c.souche)]).
            where(self.mlst.c.souche != ref)
        ).fetchall()
        return [r[0] for r in res]

    def get_all_genes(self, ref):
        res = self.connection.execute(
            select([distinct(self.mlst.c.gene)]).
            where(self.mlst.c.souche == ref)
        ).fetchall()
        return [r[0] for r in res]

    def count_sequences_per_gene(self, ref):
        res = self.connection.execute(
            select([self.mlst.c.gene, count(distinct(self.mlst.c.seqid))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_souches_per_gene(self, ref):
        res = self.connection.execute(
            select([self.mlst.c.gene, count(distinct(self.mlst.c.souche))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_genes_per_souche(self, valid_shema):
        res = self.connection.execute(
            select([self.mlst.c.souche, count(distinct(self.mlst.c.gene))])
            .where(in_(self.mlst.c.gene, valid_shema))
            .group_by(self.mlst.c.souche)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def get_sequences_number(self, ref):
        return self.connection.execute(
               select([count(distinct(self.mlst.c.seqid))])
               .where(self.mlst.c.souche != ref)
        ).fetchone()[0]

    def get_strains_distances(self, ref, valid_schema):
        m1 = self.mlst.alias()
        m2 = self.mlst.alias()

        res = self.connection.execute(
            select(
                [m1.c.souche, m2.c.souche, count(distinct(m1.c.gene))])
            .select_from(
                m1.join(m2,
                        and_(
                            m1.c.gene == m2.c.gene,
                            m1.c.souche != m2.c.souche,
                            m1.c.seqid != m2.c.seqid)))
            .where(and_(
                in_(m1.c.gene, valid_schema),
                m1.c.souche != ref,
                m2.c.souche != ref))
            .group_by(m1.c.souche, m2.c.souche)
        ).fetchall()

        distance = {}
        for r in res:
            x = distance.setdefault(r[0], {})
            x[r[1]] = r[2]

        return distance

    def get_mlst(self, ref, valid_schema):
        res = self.connection.execute(
            select([self.mlst.c.gene, self.mlst.c.souche, func.group_concat(self.mlst.c.seqid, ';')])
            .where(and_(self.mlst.c.souche != ref,
                        in_(self.mlst.c.gene, valid_schema)))
            .group_by(self.mlst.c.gene, self.mlst.c.souche)
        ).fetchall()

        mlst = {}

        for r in res:
            x = mlst.setdefault(r[0], {})
            x[r[1]] = r[2]
        return mlst

    def commit(self):
        self.transaction.commit()

    def rollback(self):
        self.transaction.rollback()

    def close(self):
        self.engine.dispose()
