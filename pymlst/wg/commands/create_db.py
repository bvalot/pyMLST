"""create_db CLI command file."""

import os
import click

from pymlst import open_wg


@click.command()
@click.option('--concatenate', '-c',
              is_flag=True,
              help='Automatically concatenate genes with duplicated sequences')
@click.option('--remove', '-r',
              is_flag=True,
              help='Automatically remove genes with duplicated sequences')
@click.argument('coregene', type=click.File('r'))
@click.argument('database', type=click.File('w'))
def cli(coregene, database, concatenate, remove):
    """Create a wgMLST database from a template"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.create(coregene, concatenate, remove)
