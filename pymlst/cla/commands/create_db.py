"""create_db CLI command file."""

import os
import click

import pymlst


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

    with pymlst.open_cla(os.path.abspath(database.name)) as mlst:
        mlst.create(scheme, alleles)
