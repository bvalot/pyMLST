"""Initial

Revision ID: 1f96d027f4aa
Revises: 
Create Date: 2024-04-29 10:11:29.815236

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '1f96d027f4aa'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    
    engine = op.get_bind()
    inspector = sa.inspect(engine)
    tables = inspector.get_table_names()

    if 'typer' not in tables:
        op.create_table('typer',
            sa.Column('id', sa.Integer(), nullable=False),
            sa.Column('sequence', sa.Text(), nullable=False),
            sa.Column('typing', sa.Text(), nullable=True),
            sa.Column('allele', sa.Text(), nullable=False),
            sa.PrimaryKeyConstraint('id'),
            sa.UniqueConstraint('sequence'))


    if 'mlst_type' not in tables:
        table = op.create_table('mlst_type',
            sa.Column('name', sa.String(7), nullable=False,
                      primary_key=True))
        data = [ { 'name' : 'pytyper'}]
        op.bulk_insert(table, data)


def downgrade():
    op.drop_table(['typer', 'mlst_type'])

