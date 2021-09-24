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
