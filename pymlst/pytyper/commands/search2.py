""" Search CLI command file for raw reads """

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command(name='search2')
@click.option('--identity', '-i',
              type=click.FLOAT,
              help='Minimum identity to search gene (default=0.9).')
@click.option('--coverage', '-c',
              type=click.FLOAT,
              help='Minimum coverage to search gene (default=0.95).')
@click.option('--reads', '-r',
              type=click.INT,
              help='Minimum reads coverage to search gene (default=10).')
@click.option('--paired/--single', default=True, 
              help= "Defines type of fastqs files.")
@click.option('--fasta', '-f',
              type=click.File('w'),
              help='Writes fasta file with gene allele.')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Writes search result to (default:stdout).')
@click.argument("method",
                type=click.Choice(['fim', 'spa', 'clmt']),
                required=True)
@click.argument('fastqs',
                type=click.File('r'),
                required=True,
                nargs=-1)


def cli(method, fastqs, **kwargs):
    """Searches strain type using specified METHOD from FASTQS(.gz) raw reads."""

    try:
        with pymlst.open_typer(method) as typer:
            typer.multi_read(fastqs, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))

