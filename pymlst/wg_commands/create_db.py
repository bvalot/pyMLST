#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a wgMLST database"""

import os
import click

from pymlst.api.core import open_wg


@click.command(name="create_db")
@click.option('--concatenate', '-c',
              is_flag=True,
              help='Automatically concatenate genes with duplicated sequences')
@click.option('--remove', '-r',
              is_flag=True,
              help='Automatically remove genes with duplicated sequences')
@click.argument('coregene', type=click.File('r'))
@click.argument('database', type=click.File('w'))
def cli(coregene, database, concatenate, remove):
    """Create a wgMLST database from a template"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.create(coregene, concatenate, remove)
