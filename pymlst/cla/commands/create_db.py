"""create_db CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions


@click.command()
@click.argument('database',
                type=click.File('w'))
@click.argument('scheme',
                type=click.File('r'))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)
def cli(database, scheme, alleles):
    """Create a classical MLST database from a sheme"""

    database.close()

    try:

        with pymlst.open_cla(os.path.abspath(database.name)) as mlst:
            mlst.create(scheme, alleles)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
