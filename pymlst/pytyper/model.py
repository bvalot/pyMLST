from sqlalchemy import MetaData, Column, Table, Integer, Text, ForeignKey, Index

metadata = MetaData()

typer = Table("typer", metadata,
			Column("id", Integer, primary_key=True),
			Column("sequence", Text, unique=True),
			Column("typing", Text),
			Column("allele", Text))


