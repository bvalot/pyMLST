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
@click.option('--list_file', '-l',
              type=click.File('r'),
              help='File containing list of coregenes to extract (default:all coregenes).')
@click.argument('database',
                type=click.File('r'))
def cli(database, **kwargs):
    """Extract sequences from a wgMLST DATABASE."""

    database.close()

    seq_kwargs, out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(SequenceExtractor(**seq_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
