"""extract strains CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import StrainExtractor, TableExtractorCommand

@click.command(name='strain', cls=TableExtractorCommand)
@click.option('--count', '-c',
              is_flag=True,
              help='Count the number of gene present in the database for each strains.')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export strain list to (default=stdout).')
@click.argument('database', type=click.Path(exists=True))
def cli(database, **kwargs):
    """Extracts a list of strains from a wgMLST DATABASE."""

    tab_kwargs,out_kwargs = utils.get_output(utils.clean_kwargs(kwargs))

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.extract(StrainExtractor(**tab_kwargs), **out_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
