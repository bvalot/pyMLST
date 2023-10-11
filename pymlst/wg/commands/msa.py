"""Multiple Sequence Alignment CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import MsaExtractor


@click.command(name='msa')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output result in fasta format (default:stdout).')
@click.option('--file', '-f',
              type=click.File('r'),
              help='file containing list of coregenes to extract (default:all coregenes).')
@click.option('--realign', '-r',
              is_flag=True,
              help='Realigns genes with same length (Default:No).')
@click.argument('database',
                type=click.Path(exists=True))
def cli(database, **kwargs):
    """Computes Multiple Sequence Alignment from a wgMLST DATABASE."""

    seq_kwargs, out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.extract(MsaExtractor(**seq_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
