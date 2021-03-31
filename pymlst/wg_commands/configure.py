"""Configure executables paths"""

import click

from pymlst import update_binary_paths
from pymlst import reset_config
from pymlst import list_binary_paths


@click.command()
@click.option('--blat', '-b',
              type=click.Path(exists=True, dir_okay=False),
              help='Blat executable absolute path')
@click.option('--mafft', '-m',
              type=click.Path(exists=True, dir_okay=False),
              help='Mafft executable absolute path')
@click.option('--reset', '-r',
              is_flag=True,
              help='Reset the configuration')
def cli(blat, mafft, reset):
    """Configure executables paths"""
    if reset:
        reset_config()
        click.echo('Resetting the configuration...')

    if mafft or blat:
        paths = {}
        if blat:
            paths['blat'] = blat
        if mafft:
            paths['mafft'] = mafft
        update_binary_paths(paths)

    config = list_binary_paths()
    if len(config) > 0:
        click.echo('--- Configuration ---')
        for k, v in config:
            click.echo(k + ': ' + v)
        click.echo('---------------------')
    elif not reset:
        click.echo('The configuration is empty')
