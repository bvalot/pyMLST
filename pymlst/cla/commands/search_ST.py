"""search_ST CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command(name='search_ST')
@click.option('--identity', '-i',
              type=click.FLOAT,
              help='Minimum identity to search gene (default=0.9).')
@click.option('--coverage', '-c',
              type=click.FLOAT,
              help='Minimum coverage to search gene (default=0.9).')
@click.option('--fasta', '-f',
              type=click.File('w'),
              help='Write fasta file with gene allele.')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Write ST search result to (default:stdout).')
@click.argument('database',
                type=click.Path(exists=True))
@click.argument('genome',
                type=click.File('r'))


def cli(genome, database, **kwargs):
    """Search ST number for an assembly GENOME using an mlst DATABASE."""


    try:

        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.search_st(genome, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
