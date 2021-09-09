"""find_subgraph CLI command file."""

import click

from pymlst.common import utils, exceptions
from pymlst.wg import core

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(name='find_subgraph',context_settings=CONTEXT_SETTINGS)
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output group files (default:stdout).')
@click.option('--threshold', '-t',
              type=click.INT,
              help='Minimum distance to conserve '
                   'for extraction of group (default:50).')
@click.option('--export', '-e',
              type=click.STRING,
              help='Export type (default:group).')
@click.argument('distance',
                type=click.File('r'))
def cli(distance, **kwargs):
    """Search group of strains at a DISTANCE threshold."""

    try:

        core.find_subgraph(distance, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
