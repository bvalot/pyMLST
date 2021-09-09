"""extract_table CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import TableExtractor, ExportType, ExtractorCommand


@click.command(name='extract_table', cls=ExtractorCommand)
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export MLST table to (default=stdout).')
@click.option('--export', '-e',
              type=click.Choice(ExportType.list_types()),
              help='Defined the export format.')
@click.option('--count', '-c',
              is_flag=True,
              help='In strain mode, count the number of gene present in the database.')
@click.argument('database', type=click.File('r'))
def cli(database, **kwargs):
    """Extract an MLST table from a wgMLST DATABASE."""

    database.close()

    tab_kwargs = utils.clean_kwargs(kwargs)

    if 'output' in tab_kwargs:
        ext_kwargs = {'output': tab_kwargs['output']}
        tab_kwargs.pop('output')
    else:
        ext_kwargs = {}

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(TableExtractor(**tab_kwargs), **ext_kwargs)

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
