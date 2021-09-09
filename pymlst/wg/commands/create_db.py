"""create_db CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(name='create_db',context_settings=CONTEXT_SETTINGS)
@click.option('--concatenate', '-c',
              is_flag=True,
              help='Automatically concatenate genes with duplicated sequences.')
@click.option('--remove', '-r',
              is_flag=True,
              help='Automatically remove genes with duplicated sequences.')
@click.argument('database', type=click.File('w'))
@click.argument('coregene', type=click.File('r'))

def cli(database, coregene, concatenate, remove):
    """Create a wgMLST DATABASE from a template COREGENE."""

    database.close()

    try:

        with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
            mlst.create(coregene, concatenate, remove)

    except exceptions.DuplicatedGeneSequence as err:
        raise click.UsageError('{}, use -c or -r options to manage it'
                               .format(str(err)))
    except exceptions.PyMLSTError as err:
        raise click.UsageError(str(err))
