"""extract_sequence CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import SequenceExtractor


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output result in fasta format (default:stdout)')
@click.option('--list', '-l',
              type=click.File('r'),
              help='List of coregenes to extract (default:all)')
@click.option('--align', '-a',
              is_flag=True,
              help='Report a concatened multi-fasta file '
              'instead of only gene files (default:No)')
@click.option('--realign', '-r',
              is_flag=True,
              help='Realign genes with same length (Default:No)')
@click.option('--mincover', '-m',
              type=click.INT,
              help='Minimun number of strain found '
              'to keep a coregene (default:1)')
@click.argument('database',
                type=click.File('r'))
def cli(database, list, **kwargs):
    """Get sequences from a wgMLST database"""

    database.close()

    seq_kwargs = utils.clean_kwargs(kwargs)

    if 'output' in seq_kwargs:
        ext_kwargs = {'output': seq_kwargs['output']}
        seq_kwargs.pop('output')
    else:
        ext_kwargs = {}

    try:
        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(SequenceExtractor(list, **seq_kwargs), **ext_kwargs)
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
