#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Remove strains to an wgMLST database"""
import os

import click

from pymlst.api.core import open_wg

from pymlst.lib.benchmark import benchmark

desc = "Remove strain to a wgMLST database and sequences specificaly associated"


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of strains to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('strains',
                type=str, nargs=-1)
@benchmark
def cli(list, database, strains):
    """Remove strain to a wgMLST database and sequences specificaly associated"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.remove_strain(strains, list)
