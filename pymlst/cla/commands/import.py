"""import CLI command file."""

import logging
import os
import sys
import tempfile

import click

import requests

import pymlst

from pymlst.common import web, exceptions
from pymlst.common import utils


@click.command(name='import')
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
        raise click.ClickException('Could not retrieve online data')
    except requests.exceptions.ConnectionError:
        raise click.ClickException('Could not access to the server, please verify your internet connection')
    except requests.exceptions.Timeout:
        raise click.ClickException('The server took too long to respond')
    except web.StructureError:
        raise click.ClickException('It seems like the structure of the website/API changed '
                                   'since this application was developed.')
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
