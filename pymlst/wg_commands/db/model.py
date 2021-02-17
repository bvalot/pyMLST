from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, Text


Base = declarative_base()


class Sequence(Base):
    __tablename__ = 'sequences'

    id = Column(Integer, primary_key=True)
    sequence = Column(Text, unique=True)


class Mlst(Base):
    __tablename__ = 'mlst'

    id = Column(Integer, primary_key=True)
    souche = Column(Text)
    gene = Column(Text)
    seqid = Column(Integer)
