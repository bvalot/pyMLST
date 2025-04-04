"""initial

Revision ID: 52ae99cb5f33
Revises: 
Create Date: 2021-05-21 10:23:49.557993

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '52ae99cb5f33'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # This is the initial revision created after the PyMLST refactoring.
    # The old un-versioned databases data are untouched.
    # Old databased indexes are dropped and replaced by new ones.
    # A new alembic_version table is added automatically to enable versioning.

    engine = op.get_bind()
    inspector = sa.inspect(engine)
    tables = inspector.get_table_names()

    if 'sequences' not in tables:
        op.create_table('sequences',
            sa.Column('id', sa.Integer(), nullable=False),
            sa.Column('sequence', sa.Text(), nullable=True),
            sa.PrimaryKeyConstraint('id'),
            sa.UniqueConstraint('sequence'))

    if 'mlst' not in tables:
        op.create_table('mlst',
            sa.Column('id', sa.Integer(), nullable=False),
            sa.Column('souche', sa.Text(), nullable=True),
            sa.Column('gene', sa.Text(), nullable=True),
            sa.Column('seqid', sa.Integer(), nullable=True),
            sa.ForeignKeyConstraint(['seqid'], ['sequences.id'], ),
            sa.PrimaryKeyConstraint('id'))

    if 'mlst_type' not in tables:
        table = op.create_table('mlst_type',
                    sa.Column('name', sa.String(length=4), nullable=False,
                              primary_key=True))
        data = [ { 'name':  'wg' } ]
        op.bulk_insert(table, data)

    indexes = inspector.get_indexes('mlst')
    for ind in indexes:
        op.drop_index(ind['name'])

    op.create_index('ix_gene', 'mlst', ['gene'], unique=False)
    op.create_index('ix_seqid', 'mlst', ['seqid'], unique=False)
    op.create_index('ix_souche', 'mlst', ['souche'], unique=False)
    op.create_index('ix_souche_gene_seqid', 'mlst', ['gene', 'souche', 'seqid'], unique=False)


def downgrade():
    # Remove index and mlst_type
    op.drop_index('ix_souche_gene_seqid', table_name='mlst')
    op.drop_index('ix_souche', table_name='mlst')
    op.drop_index('ix_seqid', table_name='mlst')
    op.drop_index('ix_gene', table_name='mlst')
    op.drop_table('mlst_type')

    # Rebuild older index
    op.create_index('ID_gene', 'mlst', ['gene'], unique=False)
    op.create_index('ID_seqid', 'mlst', ['seqid'], unique=False)
    op.create_index('ID_souche', 'mlst', ['souche'], unique=False)
    op.create_index('ID_index', 'mlst', ['souche', 'gene'], unique=False)    
