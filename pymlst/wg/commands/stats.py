"""extract stats CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import StatsExtractor

@click.command(name='stats')

@click.argument('database', type=click.File('r'))
def cli(database, **kwargs):
    """Extract stats from a wgMLST DATABASE."""

    database.close()
    
    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.extract(StatsExtractor())

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
