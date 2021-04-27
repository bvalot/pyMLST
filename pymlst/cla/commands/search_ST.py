"""search_ST CLI command file."""

import sys
import os

import click

from pymlst import open_cla


@click.command()
@click.option('--identity', '-i',
              type=float, default=0.9,
              help='Minimum identity to search gene (default=0.9)')
@click.option('--coverage', '-c',
              type=float, default=0.9,
              help='Minimum coverage to search gene (default=0.9)')
@click.option('--fasta', '-f',
              type=click.File('w'),
              help='Write fasta file with gene allele')
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Write ST search result to (default:stdout)')
@click.argument('genome',
                type=click.File('r'))
@click.argument('database',
                type=click.File('r'))
def cli(identity, coverage, fasta, output, genome, database):
    """Search ST number for an assembly"""

    database.close()

    with open_cla(os.path.abspath(database.name)) as mlst:
        mlst.search_st(genome, identity, coverage, fasta, output)
