"""recombinaison CLI command file."""

import click

from pymlst.common import utils, exceptions
from pymlst.wg import core


@click.command(name='recombinaison')
@click.option('--output', '-o',
              type=click.File('w'),
              help='Output number of variations by genes (default:stdout).')
@click.argument('genes',
                type=click.File('r'))
@click.argument('alignment',
                type=click.File('r'))
def cli(genes, alignment, **kwargs):
    """Search potential gene re-combinations from wgMLST database export."""

    try:

        core.find_recombination(genes, alignment, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
