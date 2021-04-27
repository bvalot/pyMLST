"""extract_sequence CLI command file."""

import sys
import os

import click

from pymlst import open_wg
from pymlst.wg.extractors import SequenceExtractor


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output result in fasta format (default:stdout)')
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='List of coregenes to extract (default:all)')
@click.option('--align', '-a',
              is_flag=True,
              help='Report a concatened multi-fasta file '
              'instead of only gene files (default:No)')
@click.option('--realign', '-r',
              is_flag=True,
              help='Realign genes with same length (Default:No)')
@click.option('--mincover', '-m',
              type=int, default=1,
              help='Minimun number of strain found '
              'to keep a coregene (default:1)')
@click.argument('database',
                type=click.File('r'))
def cli(output, list, align, realign, mincover, database):
    """Get sequences from a wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.extract(SequenceExtractor(list, align, realign, mincover), output)
