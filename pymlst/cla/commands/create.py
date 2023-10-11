"""create CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions


@click.command(name='create')
@click.option('--force', '-f',
              is_flag=True,
              help='Overwrites alrealdy existing DATABASE')
@click.argument('database',
                type=click.Path(exists=False))
@click.argument('scheme',
                type=click.File('r'))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)


def cli(force, database, scheme, alleles):
    """Creates a classical MLST DATABASE from a SCHEME csv and ALLELES files."""

    try:

        if os.path.exists(database):
            if force:
                open(database, "w").close()
            else:
                raise exceptions.PyMLSTError("Database alreadly exists, use --force to override it")

        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.create(scheme, alleles)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
