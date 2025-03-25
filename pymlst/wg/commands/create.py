"""create CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions, utils

@click.command(name='create')
@click.option('--force', '-f',
              is_flag=True,
              help='Overwrite alrealdy existing DATABASE')
@click.option('--concatenate', '-c',
              is_flag=True,
              help='Automatically concatenates GENES with duplicated sequences.')
@click.option('--remove', '-r',
              is_flag=True,
              help='Automatically removes GENES with duplicated sequences.')
@click.option('--species', '-s',
              type=click.STRING,
              help='Name of the species')
@click.option('--version', '-V',
              type=click.STRING,
              help='Version of the database')
@click.argument('database', type=click.Path(exists=False))
@click.argument('coregene', type=click.File('r'))

def cli(force, species, version, database, **kwargs):
    """Creates a wgMLST DATABASE from a template COREGENE."""
       
    try:

        if os.path.exists(database):
            if force:
                open(database, "w").close()
            else:
                raise exceptions.PyMLSTError("Database alreadly exists, use --force to override it")

        with pymlst.open_wg(os.path.abspath(database)) as mlst:
            mlst.create(**utils.clean_kwargs(kwargs))
            mlst.add_infos("custom", species, version)

    except exceptions.DuplicatedGeneSequence as err:
        raise click.UsageError('{}, use -c or -r options to manage it'
                               .format(str(err)))
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
