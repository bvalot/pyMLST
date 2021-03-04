from sqlalchemy import create_engine, MetaData, Table, Column, Integer, Text


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

    def commit(self):
        self.transaction.commit()

    def rollback(self):
        self.transaction.rollback()

    def close(self):
        self.engine.dispose()