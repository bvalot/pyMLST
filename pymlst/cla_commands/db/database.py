from sqlalchemy import create_engine, MetaData, Table, Column, Integer, Text, select, distinct, and_


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

    def add_sequence_safe(self, sequence):
        """Inserts a sequence if it doesn't already exist yet"""
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