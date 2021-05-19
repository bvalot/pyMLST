"""remove_strain CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command()
@click.option('--list', '-l',
              type=click.File('r'),
              help='File list of strains to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('strains',
                type=str, nargs=-1)
def cli(database, strains, **kwargs):
    """Remove strain to a wgMLST database and sequences specificaly associated"""

    database.close()

    try:
        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.remove_strain(strains, **utils.clean_kwargs(kwargs))
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
