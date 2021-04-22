#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL
import logging
import os
import click
import tempfile

import requests

from pymlst import open_wg
from pymlst.common import utils
from pymlst.common.web import prompt_cgmlst, build_coregene


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

    utils.create_logger()

    try:

        url = prompt_cgmlst(' '.join(species), prompt)
        if url == '':
            logging.info('No choice selected')
            return
        logging.info('Downloading the core genome...')

        with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:

            skipped = build_coregene(url, tmp)
            tmp.close()
            if len(skipped) > 0:
                logging.info('Skipped the following malformed file(s): ' + ', '.join(skipped))

            with open_wg(os.path.abspath(database.name)) as mlst:
                mlst.create(tmp.name)

    except requests.exceptions.HTTPError:
        logging.error('An error occurred while retrieving online data')
        exit(1)
    except requests.exceptions.ConnectionError:
        logging.error('Couldn\'t access to the server, please verify your internet connection')
        exit(2)
    except requests.exceptions.Timeout:
        logging.error('The server took too long to respond')
        exit(3)
