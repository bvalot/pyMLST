"""extract MLST table CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import MlstExtractor, TableExtractorCommand

@click.command(name='mlst', cls=TableExtractorCommand)
@click.option('--form', '-f',
              type=click.Choice(["default", "grapetree"]),
              help='Specify format of output')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export strain list to (default=stdout).')
@click.argument('database', type=click.File('r'))
def cli(database, **kwargs):
    """Extract an MLST table from a wgMLST DATABASE."""

    database.close()

    tab_kwargs,out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))
    
    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(MlstExtractor(**tab_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
