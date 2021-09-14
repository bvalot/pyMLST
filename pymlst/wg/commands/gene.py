"""extract gene CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import GeneExtractor, TableExtractorCommand

@click.command(name='gene', cls=TableExtractorCommand)
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export MLST table to (default=stdout).')
@click.argument('database', type=click.File('r'))
def cli(database, **kwargs):
    """Extract an genes list from a wgMLST DATABASE."""

    database.close()

    tab_kwargs,out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(GeneExtractor(**tab_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
