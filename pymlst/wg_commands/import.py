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

import requests

from pymlst.api.core import open_wg, create_logger
from pymlst.lib.web import prompt_cgmlst, build_coregene


@click.command()
@click.option('--prompt/--no-prompt',
              default=True)
@click.argument('database',
                type=click.File('w'))
@click.argument('species',
                type=click.STRING,
                nargs=-1)
def cli(prompt, database, species):
    """Create a wgMLST database from an online resource"""

    logger = create_logger()

    try:

        url = prompt_cgmlst(' '.join(species), prompt)
        if url == '':
            logger.info('No choice selected')
            return
        logger.info('Downloading the core genome...')

        with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:

            skipped = build_coregene(url, tmp)
            tmp.close()
            if len(skipped) > 0:
                logger.info('Skipped the following malformed file(s): ' + ', '.join(skipped))

            with open_wg(os.path.abspath(database.name)) as mlst:
                mlst.create(tmp.name)

    except requests.exceptions.HTTPError:
        logger.error('An error occurred while retrieving online data')
        exit(1)
    except requests.exceptions.ConnectionError:
        logger.error('Couldn\'t access to the server, please verify your internet connection')
        exit(2)
    except requests.exceptions.Timeout:
        logger.error('The server took too long to respond')
        exit(3)
