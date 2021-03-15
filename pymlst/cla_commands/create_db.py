#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a classical MLST database"""
import os
import click

from pymlst.api.core import open_cla

from pymlst.lib.benchmark import benchmark

desc = "Create a classical MLST database from a sheme"


@click.command()
@click.argument('database',
                type=click.File('w'))
@click.argument('scheme',
                type=click.File('r'))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)
@benchmark
def cli(database, scheme, alleles):
    """Create a classical MLST database from a sheme"""

    database.close()
    
    with open_cla(os.path.abspath(database.name)) as mlst:
        mlst.create(scheme, alleles)
