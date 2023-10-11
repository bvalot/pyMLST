"""extract stats CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions
from pymlst.wg.extractors import StatsExtractor

@click.command(name='stats')

@click.argument('database', type=click.Path(exists=True))
def cli(database, **kwargs):
    """Extracts stats from a wgMLST DATABASE."""

    try:

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.extract(StatsExtractor())

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
