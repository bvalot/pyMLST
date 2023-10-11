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
@click.option('--force', '-f',
              is_flag=True,
              help='Overwrites alrealdy existing DATABASE')
@click.option('--prompt/--no-prompt',
              default=True,
              help='Do not prompt if multiple '
                   'choices are found, fail instead.')
@click.option('--mlst', '-m',
              type=click.STRING,
              default='',
              help='Specifies the desired MLST scheme name.')
@click.argument('database',
                type=click.Path(exists=False))
@click.argument('species',
                type=click.STRING,
                nargs=-1)


def cli(force, prompt, mlst, database, species):
    """Creates a claMLST DATABASE from an online resource.

    The research can be filtered by adding a SPECIES name."""

    utils.create_logger()

    try:

        if os.path.exists(database):
            if force:
                open(database, "w").close()
            else:
                raise exceptions.PyMLSTError("Database alreadly exists, use --force to override it")

        url = web.retrieve_mlst(' '.join(species), prompt, mlst)

        if url is None:
            logging.info('No choice selected')
            return

        logging.info('Downloading mlst...')

        with tempfile.TemporaryDirectory() as tmp_dir, \
                pymlst.open_cla(os.path.abspath(database)) as mlst_db:

            web.get_mlst_files(url, tmp_dir)

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
