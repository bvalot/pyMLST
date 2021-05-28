from sqlalchemy import MetaData, Column, Table, Integer, Text, ForeignKey, Index

metadata = MetaData()

sequences = Table('sequences', metadata,
                  Column('id', Integer, primary_key=True),
                  Column('sequence', Text, unique=True))

mlst = Table('mlst', metadata,
             Column('id', Integer, primary_key=True),
             Column('souche', Text),
             Column('gene', Text),
             Column('seqid', Integer, ForeignKey(sequences.c.id)),
             Index('ix_souche', 'souche'),
             Index('ix_gene', 'gene'),
             Index('ix_seqid', 'seqid'),
             Index('ix_souche_gene_seqid', 'gene', 'souche', 'seqid'))
