"""add_strain CLI command file."""

import os

import click

from pymlst import open_wg


@click.command()
@click.option('--strain', '-s',
              type=str, default=None,
              help='Name of the strain (default:genome name)')
@click.option('--identity', '-i',
              type=float, default=0.95,
              help='Minimum identity to search gene (default=0.95)')
@click.option('--coverage', '-c',
              type=float, default=0.9,
              help='Minimum coverage to search gene (default=0.9)')
@click.argument('genome',
                type=click.File("r"))
@click.argument('database',
                type=click.File("r"))
def cli(strain, identity, coverage, genome, database):
    """Add a strain to the wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.add_strain(genome, strain, identity, coverage)
