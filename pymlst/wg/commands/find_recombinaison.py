"""find_recombinaison CLI command file."""

import click

from pymlst.common import utils, exceptions
from pymlst.wg import core

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(name='find_recombinaison',context_settings=CONTEXT_SETTINGS)
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
