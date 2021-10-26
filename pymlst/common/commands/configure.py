"""configure CLI command file."""

import click

from pymlst import config

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--blat', '-b',
              type=click.Path(exists=True, dir_okay=False),
              help='Blat executable absolute path.')
@click.option('--kma', '-k',
              type=click.Path(exists=True, dir_okay=False),
              help='Kma executable absolute path.')
@click.option('--mafft', '-m',
              type=click.Path(exists=True, dir_okay=False),
              help='Mafft executable absolute path.')
@click.option('--log', '-l',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='Level of logging, default=INFO')
@click.option('--reset', '-r',
              is_flag=True,
              help='Reset the configuration.')
def cli(blat, kma, mafft, log, reset):
    """Configure executables paths and log level."""
    if reset:
        config.reset_binary_paths()
        config.set_logging_level("INFO")
        click.echo('Resetting the configuration...')

    if mafft or blat or kma:
        paths = {}
        if blat:
            paths['blat'] = blat
        if kma:
            paths['kma'] = kma
        if mafft:
            paths['mafft'] = mafft
        config.update_binary_paths(paths)
    if log:
        config.set_logging_level(log)

    paths = config.list_binary_paths()
    log = config.get_logging_level()
    click.echo('--- Configuration ---')
    if len(paths) > 0:
        for key, value in paths:
            click.echo(key + ': ' + value)
        click.echo('---------------------')
    click.echo('LOG : ' + log)
    click.echo('---------------------')
