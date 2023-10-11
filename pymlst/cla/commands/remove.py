"""remove CLI command file."""

import os

import click

import pymlst
from pymlst.common import utils, exceptions


@click.command(name='remove')
@click.argument('database',
                type=click.Path(exists=True))
@click.argument('gene',
                type=click.STRING)
@click.argument('allele',
                type=click.INT)


def cli(database, **kwargs):
    """Removes ALLELE sequence from the GENE on a mlst DATABASE."""
    
    try:
        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.remove_allele(**utils.clean_kwargs(kwargs))
                
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
