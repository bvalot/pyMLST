""" search CLI command file. """

import os
import click

import pymlst
from pymlst.pytyper import model
from pymlst.pytyper.method import FIM, SPA, CLMT
from pymlst.common import utils, exceptions

@click.command(name="search")
@click.option("--identity", "-i",
              type=click.FLOAT,
              help="Minimum identity to search gene (default=0.9).")
@click.option("--coverage", "-c",
              type=click.FLOAT,
              help="Minimum coverage to search gene (default=0.9).")
@click.option('--fasta', '-f',
              type=click.File('w'),
              help='Writes fasta file with gene allele.')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Writes search result to (default:stdout).')

# Database is initialized automatically without intervention from user
@click.argument('method',
                type=click.Choice([FIM, SPA, CLMT]),
                required=True)
@click.argument('genomes',
                type=click.File('r'),
                required=True,
                nargs=-1)

def cli(method, genomes, **kwargs):
    """Searches strain type using specified METHOD for an assembly GENOME.
    
    ...

    """
    
    try:
        with pymlst.open_typer(method) as typer:
            typer.multi_search(genomes, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))


