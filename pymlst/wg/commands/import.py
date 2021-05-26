"""import CLI command file."""

import logging
import os
import tempfile

import click
import requests

import pymlst
from pymlst.common import utils, web, exceptions


@click.command(name='import')
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
