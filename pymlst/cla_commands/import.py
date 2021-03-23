#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Initialize a database from an online base"""
import os
import click
import tempfile

from pymlst.api.core import open_cla

from pymlst.lib.benchmark import benchmark
from pymlst.lib.web import get_mlst_files
from pymlst.api.core import create_logger

desc = "Initialize a database from an online base"


@click.command()
@click.argument('species',
                type=str)
@click.argument('database',
                type=click.File('w'))
def cli(species, database):
    """Initialize a database from an online base"""

    database.close()

    with tempfile.TemporaryDirectory() as tmp_dir, \
            open_cla(os.path.abspath(database.name)) as mlst:

        create_logger().info('Downloading the mlst...')
        get_mlst_files(species, tmp_dir)

        mlst.create(open(tmp_dir + '/profiles.csv', 'rt'),
                    [open(tmp_dir + '/locus/' + locus, 'r') for locus in os.listdir(tmp_dir + '/locus')])
