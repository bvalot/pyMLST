
import os
import click

import pymlst
from pymlst.common import utils, exceptions


@click.command(name="remove")
@click.option('--item', '-i', default='strains', show_default=True,
              type=click.Choice(['strains','genes'], case_sensitive=False),
              help= "Choose the item you wish to remove : strain or genes")

@click.option('--file', '-f',type=click.File('r'),
              help='File list of genes or strains to removed on the wgMLST database.')

@click.argument('database', type=click.File('r'), nargs=1)
@click.argument('genes_or_strains', required=False, type=str, nargs=-1)


def cli(database, item, genes_or_strains, **kwargs):
    """Remove STRAINS or GENES from a wgMLST DATABASE."""
    
    database.close()

    if item == "strains":
        click.echo("We will remove one or more strain(s)")
        try:
            with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
                mlst.remove_strain(genes_or_strains, **utils.clean_kwargs(kwargs))
                
        except exceptions.PyMLSTError as err:
            raise click.ClickException(str(err))
       
    if item == "genes":
        click.echo("We will remove one or more gene(s)")
        try:
            with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
                mlst.remove_gene(genes_or_strains, **utils.clean_kwargs(kwargs))

        except exceptions.PyMLSTError as err:
            raise click.ClickException(str(err))