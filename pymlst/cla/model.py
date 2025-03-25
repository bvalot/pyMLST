from sqlalchemy import MetaData, Table, Column, Integer, Text

metadata = MetaData()

sequences = Table('sequences', metadata,
                  Column('id', Integer, primary_key=True),
                  Column('sequence', Text, unique=True),
                  Column('gene', Text),
                  Column('allele', Integer))

mlst = Table('mlst', metadata,
             Column('id', Integer, primary_key=True),
             Column('st', Integer),
             Column('gene', Text),
             Column('allele', Integer))

mlst_type = Table('mlst_type', metadata,
                  Column('name', Text),
                  Column('source', Text),
                  Column('species', Text),
                  Column('mlst', Text),
                  Column('version', Text))
