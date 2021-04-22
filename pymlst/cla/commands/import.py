#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Initialize a database from an online base"""
import logging
import os
import click
import tempfile

import requests

from pymlst import open_cla

from pymlst.common.web import get_mlst_files, prompt_mlst
from pymlst.common import utils


@click.command()
@click.option('--prompt/--no-prompt',
              default=True)
@click.option('--mlst', '-m',
              type=click.STRING)
@click.argument('database',
                type=click.File('w'))
@click.argument('species',
                type=click.STRING,
                nargs=-1)
def cli(prompt, mlst, database, species):
    """Initialize a database from an online base"""

    database.close()
    utils.create_logger()

    try:

        url = prompt_mlst(' '.join(species), prompt, mlst)

        if url is None:
            logging.info('No choice selected')
            return
        else:
            logging.info('Downloading mlst...')

        with tempfile.TemporaryDirectory() as tmp_dir, \
                open_cla(os.path.abspath(database.name)) as mlst:

            get_mlst_files(tmp_dir, url=url)

            mlst.create(open(tmp_dir + '/profiles.csv', 'rt'),
                        [open(tmp_dir + '/locus/' + locus, 'r') for locus in os.listdir(tmp_dir + '/locus')])

    except requests.exceptions.HTTPError:
        logging.error('An error occurred while retrieving online data')
        exit(1)
    except requests.exceptions.ConnectionError:
        logging.error('Couldn\'t access to the server, please verify your internet connection')
        exit(2)
    except requests.exceptions.Timeout:
        logging.error('The server took too long to respond')
        exit(3)
