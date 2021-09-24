from sqlalchemy import MetaData, Table, String, Column

metadata = MetaData()

mlst_type = Table('mlst_type', metadata,
                  Column('name', String(length=4), primary_key=True))
