"""find_subgraph CLI command file."""

import click

from pymlst.common import utils
from pymlst.wg import core


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output group files (default:stdout)')
@click.option('--threshold', '-t',
              type=click.INT,
              help='Minimum distance to conserve '
                   'for extraction of group (default:50)')
@click.option('--export', '-e',
              type=click.STRING,
              help='Export type (default:group)')
@click.argument('distance',
                type=click.File('r'))
def cli(distance, **kwargs):
    """Search group os strain at a distance threshold."""

    core.find_subgraph(distance, **utils.clean_kwargs(kwargs))
