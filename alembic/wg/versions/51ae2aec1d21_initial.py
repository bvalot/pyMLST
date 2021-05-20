"""initial

Revision ID: 51ae2aec1d21
Revises: 
Create Date: 2021-05-20 15:24:56.147144

"""
from alembic import op
import sqlalchemy as sa

from sqlalchemy.engine.reflection import Inspector


# revision identifiers, used by Alembic.
revision = '51ae2aec1d21'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    engine = op.get_bind()

    metadata = sa.MetaData()

    sequences = sa.Table('sequences', metadata,
                         sa.Column('id', sa.Integer, primary_key=True),
                         sa.Column('sequence', sa.Text, unique=True))

    mlst = sa.Table('mlst', metadata,
                    sa.Column('id', sa.Integer, primary_key=True),
                    sa.Column('souche', sa.Text),
                    sa.Column('gene', sa.Text),
                    sa.Column('seqid', sa.Integer, sa.ForeignKey(sequences.c.id)))

    metadata.create_all(engine)

    sa.Index('ix_souche',
             mlst.c.souche)\
        .create(bind=engine, checkfirst=True)
    sa.Index('ix_gene',
             mlst.c.gene)\
        .create(bind=engine, checkfirst=True)
    sa.Index('ix_seqid',
             mlst.c.seqid)\
        .create(bind=engine, checkfirst=True)
    sa.Index('ix_souche_gene_seqid',
             mlst.c.gene,
             mlst.c.souche,
             mlst.c.seqid)\
        .create(bind=engine, checkfirst=True)


def downgrade():
    pass
