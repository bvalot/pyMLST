"""Info CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions, utils


@click.command(name='info')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Writes ST search result to (default:stdout).')
@click.argument('database',
                type=click.Path(exists=False))


def cli(database, **kwargs):
    """Output the information about  a classical MLST DATABASE"""

    try:
        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.get_infos(**utils.clean_kwargs(kwargs))
            
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
