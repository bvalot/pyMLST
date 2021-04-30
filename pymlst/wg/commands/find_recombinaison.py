"""find_recombinaison CLI command file."""

import click

from pymlst.common import utils
from pymlst.wg.core import find_recombination


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

    find_recombination(genes, alignment, **utils.clean_kwargs(kwargs))
