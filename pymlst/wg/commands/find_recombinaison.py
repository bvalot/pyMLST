"""find_recombinaison CLI command file."""

import click

from pymlst.common import utils
from pymlst.wg import core


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output number of variations by genes (default:stdout)')
@click.argument('genes',
                type=click.File('r'))
@click.argument('alignment',
                type=click.File('r'))
def cli(genes, alignment, **kwargs):
    """Search potential gene recombinaison from wgMLST database export"""

    core.find_recombination(genes, alignment, **utils.clean_kwargs(kwargs))
