"""Add database infos

Revision ID: a793f8f3fd83
Revises: 52ae99cb5f33
Create Date: 2025-03-14 15:38:05.090257

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'a793f8f3fd83'
down_revision = '52ae99cb5f33'
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('mlst_type', sa.Column('source', sa.String(), nullable=True)) 
    op.add_column('mlst_type', sa.Column('species', sa.String(), nullable=True))
    op.add_column('mlst_type', sa.Column('version', sa.String(), nullable=True))
    

def downgrade():
    op.drop_column('mlst_type', 'source')
    op.drop_column('mlst_type', 'species')
    op.drop_column('mlst_type', 'version')
