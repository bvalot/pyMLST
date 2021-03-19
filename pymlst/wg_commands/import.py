#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a wgMLST database from an online resource"""

import os
import click
import tempfile

from pymlst.api.core import open_wg, create_logger
from pymlst.lib.web import prompt_cgmlst, build_coregene


@click.command()
@click.argument('database', type=click.File('w'))
def cli(database):
    """Create a wgMLST database from an online resource"""

    url = prompt_cgmlst()
    create_logger().info('Downloading the core genome...')

    with tempfile.NamedTemporaryFile('w+') as tmp:

        build_coregene(url, tmp)

        with open_wg(os.path.abspath(database.name)) as mlst:
            mlst.create(tmp.name)
