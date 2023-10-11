"""add CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions

@click.command(name='add2')
@click.option('--strain', '-s',
              type=click.STRING,
              help='Name of the strain (default:genome name).')
@click.option('--identity', '-i',
              type=click.FLOAT,
              help='Minimum identity to search gene (default=0.95).')
@click.option('--coverage', '-c',
              type=click.FLOAT,
              help='Minimum coverage to search gene (default=0.9).')
@click.option('--reads', '-r',
              type=click.INT,
              help='Minimum reads coverage to search a gene (default=10).')
@click.argument('database', nargs=1, 
                type=click.Path(exists=True))
@click.argument('fastqs', nargs=-1, 
                type=click.File("r"))

def cli(fastqs, database, **kwargs):
    """Adds a strain from FASTQS(.gz) reads to the wgMLST DATABASE."""

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.add_reads(fastqs, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
