"""extract sequence CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import SequenceExtractor


@click.command(name='sequence')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output result in fasta format (default:stdout).')
@click.option('--file', '-f',
              type=click.File('r'),
              help='File containing list of coregenes to extract (default:all coregenes).')
@click.option('--reference',
              is_flag=True,
              help='Return sequence of the reference instead of strains alleles')
@click.argument('database',
                type=click.Path(exists=True))
def cli(database, **kwargs):
    """Extracts sequences from a wgMLST DATABASE."""

    seq_kwargs, out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.extract(SequenceExtractor(**seq_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
