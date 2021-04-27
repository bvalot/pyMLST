"""remove_strain CLI command file."""

import os

import click

from pymlst import open_wg


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of strains to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('strains',
                type=str, nargs=-1)
def cli(list, database, strains):
    """Remove strain to a wgMLST database and sequences specificaly associated"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.remove_strain(strains, list)
