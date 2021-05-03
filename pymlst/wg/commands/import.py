"""import CLI command file."""

import logging
import os
import sys

import click
import tempfile
import requests

import pymlst
from pymlst.common import utils
from pymlst.common import web


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

        url = web.retrieve_cgmlst(' '.join(species), prompt)
        if url == '':
            logging.info('No choice selected')
            return
        logging.info('Downloading the core genome...')

        with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:

            skipped = web.get_coregene_file(url, tmp)
            tmp.close()
            if len(skipped) > 0:
                logging.info('Skipped the following malformed file(s): %s', ', '.join(skipped))

            with pymlst.open_wg(os.path.abspath(database.name)) as mlst:
                mlst.create(tmp.name)

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
