from sqlalchemy import MetaData, Column, Table, Integer, Text, ForeignKey, Index

metadata = MetaData()

typerSeq = Table("typerSeq", metadata,
	      Column("id", Integer, primary_key=True),
	      Column("sequence", Text, unique=True),
	      Column("typing", Text),
	      Column("allele", Text))

typerSt = Table("typerSt", metadata,
	      Column("id", Integer, primary_key=True),
	      Column("typing", Text),
	      Column("st", Text),              
	      Column("allele", Text))


