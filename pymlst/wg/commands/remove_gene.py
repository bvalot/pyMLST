"""remove_gene CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command()
@click.option('--list', '-l',
              type=click.File('r'),
              help='File list of genes to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('genes',
                type=str, nargs=-1)
def cli(genes, database, **kwargs):
    """Remove gene to a wgMLST database"""

    database.close()

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.remove_gene(genes, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
