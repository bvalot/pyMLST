"""create_db CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions


@click.command(name='create_db')

@click.argument('database',
                type=click.File('w'))
@click.argument('scheme',
                type=click.Path(exists=True))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)


def cli(database, scheme, alleles):
    """Create a classical MLST DATABASE from a SCHEME csv and ALLELES files."""


    try:

        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.create(scheme, alleles)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
