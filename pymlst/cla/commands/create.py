"""create CLI command file."""

import os
import click

import pymlst
from pymlst.common import exceptions


@click.command(name='create')
@click.option('--force', '-f',
              is_flag=True,
              help='Overwrites alrealdy existing DATABASE')
@click.option('--species', '-s',
              type=click.STRING,
              help='Name of the species')
@click.option('--version', '-V',
              type=click.STRING,
              help='Version of the database')
@click.argument('database',
                type=click.Path(exists=False))
@click.argument('profile',
                type=click.File('r'))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)


def cli(force, species, version, database, profile, alleles):
    """Creates a classical MLST DATABASE from a txt PROFILE and fasta ALLELES files."""

    try:

        if os.path.exists(database):
            if force:
                open(database, "w").close()
            else:
                raise exceptions.PyMLSTError("Database alreadly exists, use --force to override it")

        with pymlst.open_cla(os.path.abspath(database)) as mlst:
            mlst.create(profile, alleles)
            mlst.add_infos("Custom", species, "", version)
            
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
