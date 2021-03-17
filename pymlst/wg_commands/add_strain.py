#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Add a strain to the wgMLST database"""

import os

import click

from pymlst.api.core import open_wg
from pymlst.lib.benchmark import benchmark

blat_path = '/usr/bin/'

desc = "Add a strain to the wgMLST database"


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
@benchmark
def cli(strain, identity, coverage, genome, database):
    """Add a strain to the wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.add_strain(genome, strain, identity, coverage)
