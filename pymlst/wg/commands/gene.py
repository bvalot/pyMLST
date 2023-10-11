"""extract gene CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import GeneExtractor, TableExtractorCommand

@click.command(name='gene', cls=TableExtractorCommand)
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export GENE list to (default=stdout).')
@click.argument('database', type=click.Path(exists=True))
def cli(database, **kwargs):
    """Extracts a list of genes from a wgMLST DATABASE."""

    tab_kwargs,out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.extract(GeneExtractor(**tab_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
