"""add CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions

@click.command(name='add')
@click.option('--strain', '-s',
              type=click.STRING,
              help='Name of the strain (default:genome name).')
@click.option('--identity', '-i',
              type=click.FLOAT,
              help='Minimum identity to search gene (default=0.95).')
@click.option('--coverage', '-c',
              type=click.FLOAT,
              help='Minimum coverage to search gene (default=0.9).')
@click.argument('database',
                type=click.Path(exists=True))
@click.argument('genome',
                type=click.File("r"))

def cli(genome, database, **kwargs):
    """Adds a strain GENOME to the wgMLST DATABASE."""

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.add_strain(genome, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
