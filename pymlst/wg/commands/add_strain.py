"""add_strain CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command(name='add_strain')
@click.option('--strain', '-s',
              type=click.STRING,
              help='Name of the strain (default:genome name).')
@click.option('--identity', '-i',
              type=click.FLOAT,
              help='Minimum identity to search gene (default=0.95).')
@click.option('--coverage', '-c',
              type=click.FLOAT,
              help='Minimum coverage to search gene (default=0.9).')
@click.argument('genome',
                type=click.File("r"))
@click.argument('database',
                type=click.File("r"))
def cli(genome, database, **kwargs):
    """Add a strain GENOME to the wgMLST DATABASE."""

    database.close()

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.add_strain(genome, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
