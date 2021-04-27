"""find_recombinaison CLI command file."""

import sys
import click

from pymlst.wg.core import find_recombination


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output number of variations by genes (default:stdout)')
@click.argument('genes',
                type=click.File('r'))
@click.argument('alignment',
                type=click.File('r'))
def cli(output, genes, alignment):
    """Search potential gene recombinaison from wgMLST database export"""

    find_recombination(genes, alignment, output)