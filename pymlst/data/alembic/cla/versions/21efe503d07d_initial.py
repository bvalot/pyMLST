"""initial

Revision ID: 21efe503d07d
Revises: 
Create Date: 2021-05-21 15:55:22.181990

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '21efe503d07d'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # This is the initial revision created after the PyMLST refactoring.
    # The old un-versioned databases data are untouched.
    # A new alembic_version table is added automatically to enable versioning.

    engine = op.get_bind()
    inspector = sa.inspect(engine)
    tables = inspector.get_table_names()

    if 'mlst' not in tables:
        op.create_table('mlst',
            sa.Column('id', sa.Integer(), nullable=False),
            sa.Column('st', sa.Integer(), nullable=True),
            sa.Column('gene', sa.Text(), nullable=True),
            sa.Column('allele', sa.Integer(), nullable=True),
            sa.PrimaryKeyConstraint('id'))

    if 'sequences' not in tables:
        op.create_table('sequences',
            sa.Column('id', sa.Integer(), nullable=False),
            sa.Column('sequence', sa.Text(), nullable=True),
            sa.Column('gene', sa.Text(), nullable=True),
            sa.Column('allele', sa.Integer(), nullable=True),
            sa.PrimaryKeyConstraint('id'),
            sa.UniqueConstraint('sequence'),
            sa.UniqueConstraint('sequence'))


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('sequences')
    op.drop_table('mlst')
    # ### end Alembic commands ###