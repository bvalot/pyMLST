"""remove_gene CLI command file."""

import os

import click

from pymlst import open_wg


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of genes to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('genes',
                type=str, nargs=-1)
def cli(list, genes, database):
    """Remove gene to a wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.remove_gene(genes, list)
