"""import CLI command file."""

import logging
import os
import sys
import tempfile

import click

import requests

import pymlst

from pymlst.common import web
from pymlst.common import utils


@click.command()
@click.option('--prompt/--no-prompt',
              default=True)
@click.option('--mlst', '-m',
              type=click.STRING,
              default='')
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

        url = web.retrieve_mlst(' '.join(species), prompt, mlst)

        if url is None:
            logging.info('No choice selected')
            return
        logging.info('Downloading mlst...')

        with tempfile.TemporaryDirectory() as tmp_dir, \
                pymlst.open_cla(os.path.abspath(database.name)) as mlst_db:

            web.get_mlst_files(tmp_dir, url=url)

            mlst_db.create(open(tmp_dir + '/profiles.csv', 'rt'),
                           [open(tmp_dir + '/locus/' + locus, 'r')
                           for locus in os.listdir(tmp_dir + '/locus')])

    except requests.exceptions.HTTPError:
        logging.error('An error occurred while retrieving online data')
        sys.exit(1)
    except requests.exceptions.ConnectionError:
        logging.error('Couldn\'t access to the server, please verify your internet connection')
        sys.exit(2)
    except requests.exceptions.Timeout:
        logging.error('The server took too long to respond')
        sys.exit(3)
    except web.StructureError:
        logging.error('It seems like the structure of the website/API changed '
                      'since this application was developed.')
        sys.exit(4)
