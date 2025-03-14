"""Add database infos

Revision ID: c0f871a99d96
Revises: 21efe503d07d
Create Date: 2025-03-14 09:29:25.322104

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'c0f871a99d96'
down_revision = '21efe503d07d'
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('mlst_type', sa.Column('source', sa.String(), nullable=True)) 
    op.add_column('mlst_type', sa.Column('species', sa.String(), nullable=True))
    op.add_column('mlst_type', sa.Column('mlst', sa.String(), nullable=True)) 
    op.add_column('mlst_type', sa.Column('version', sa.String(), nullable=True))
    


def downgrade():
    op.drop_column('mlst_type', 'source')
    op.drop_column('mlst_type', 'species')
    op.drop_column('mlst_type', 'mlst')    
    op.drop_column('mlst_type', 'version')
