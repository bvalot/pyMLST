
import os
import click
import logging

import pymlst
from pymlst.common import utils, exceptions


@click.command(name="remove")
# @click.option('--item', '-i', default='strains', show_default=True,
#               type=click.Choice(['strains','genes'], case_sensitive=False),
#               help= "Choose the item you wish to remove : strain or genes")

@click.option('--strains/--genes',
              default=True, show_default="strains", 
              help= "Choose the item you wish to remove")
@click.option('--file', '-f',type=click.File('r'),
              help='File list of genes or strains to removed on the wgMLST database.')
@click.argument('database', type=click.Path(exists=True), nargs=1)
@click.argument('genes_or_strains', required=False, type=str, nargs=-1)


def cli(database, strains, genes_or_strains, **kwargs):
    """Removes STRAINS or GENES from a wgMLST DATABASE."""

    utils.create_logger()

    try:
        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            if strains:
                logging.info("We will remove one or more strain(s)")
                mlst.remove_strain(genes_or_strains, **utils.clean_kwargs(kwargs))
       
            else :
                logging.info("We will remove one or more gene(s)")
                mlst.remove_gene(genes_or_strains, **utils.clean_kwargs(kwargs))

    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
