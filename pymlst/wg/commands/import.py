"""import CLI command file."""

import logging
import os
import tempfile

import click
import requests

import pymlst
from pymlst.common import utils, web, exceptions


@click.command(name='import')
@click.option('--force', '-f',
              is_flag=True,
              help='Overwrite alrealdy existing DATABASE')
@click.option('--prompt/--no-prompt',
              default=True,
              help='Do not prompt if multiple '
                   'choices are found, fail instead.')
@click.argument('database',
                type=click.Path(exists=False))
@click.argument('species',
                type=click.STRING,
                nargs=-1)
def cli(force, prompt, database, species):
    """Creates a wgMLST DATABASE from an online resource.

    The research can be filtered by adding a SPECIES name."""

    utils.create_logger()

    try:

        if os.path.exists(database):
            if force:
                open(database, "w").close()
            else:
                raise exceptions.PyMLSTError("Database alreadly exists, use --force to override it")
        
        url = web.retrieve_cgmlst(' '.join(species), prompt)

        if url is None:
            logging.info('No choice selected')
            return

        logging.info('Downloading the core genome...')

        with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:

            skipped = web.get_cgmlst_file(url, tmp)
            tmp.close()
            if len(skipped) > 0:
                logging.info('Skipped the following malformed file(s): %s', ', '.join(skipped))

            with pymlst.open_wg(os.path.abspath(database)) as mlst:
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
