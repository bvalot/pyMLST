"""extract_table CLI command file."""

import os
import sys

import click

from pymlst.wg.extractors import TableExtractor, ExportType
from pymlst.common import utils
from pymlst import open_wg


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              help='Export MLST table to (default=stdout)')
@click.option('--export', '-e',
              type=click.Choice(ExportType.list_types()),
              help='Defined the export format')
@click.option('--count', '-c',
              is_flag=True,
              help='In strain mode, count the number of gene present in the database')
@click.option('--mincover', '-m',
              type=click.INT,
              help='Minimun number of strain found to keep a gene (default:0)')
@click.option('--keep', '-k',
              is_flag=True,
              help='Keep only gene with different allele (omit missing)')
@click.option('--duplicate', '-d',
              is_flag=True,
              help='Conserve duplicate gene (default remove)')
@click.option('--inverse', '-V',
              is_flag=True,
              help='Keep only gene that do not ' \
              'meet the filter of mincover or keep options')
@click.argument('database', type=click.File('r'))
def cli(database, **kwargs):
    """Extract MLST table from a wgMLST database"""

    database.close()

    tab_kwargs = utils.clean_kwargs(kwargs)

    if 'output' in tab_kwargs:
        ext_kwargs = {'output': tab_kwargs['output']}
        tab_kwargs.pop('output')
    else:
        ext_kwargs = {}

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.extract(TableExtractor(**tab_kwargs), **ext_kwargs)
