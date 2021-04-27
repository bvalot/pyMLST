"""extract_table CLI command file."""

"""Extract MLST table from an wgMLST database"""
import os
import sys

import click

from pymlst.wg.extractors import TableExtractor, ExportType
from pymlst import open_wg


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Export MLST table to (default=stdout)')
@click.option('--export', '-e',
              type=click.Choice(ExportType.list_types()),
              default='mlst',
              help='Defined the export format')
@click.option('--count', '-c',
              is_flag=True,
              help='In strain mode, count the number of gene present in the database')
@click.option('--mincover', '-m',
              type=int, default=0,
              help='Minimun number of strain found to keep a gene (default:0)')
@click.option('--keep', '-k',
              is_flag=True,
              help='Keep only gene with different allele (omit missing)')
@click.option('--duplicate', '-d',
              is_flag=True, default=True,
              help='Conserve duplicate gene (default remove)')
@click.option('--inverse', '-V',
              is_flag=True,
              help='Keep only gene that do not ' \
              'meet the filter of mincover or keep options')
@click.argument('database', type=click.File('r'))
def cli(output, export, count, mincover, keep, duplicate, inverse, database):
    """Extract MLST table from a wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.extract(TableExtractor(export, count, mincover, keep, duplicate, inverse), output)
