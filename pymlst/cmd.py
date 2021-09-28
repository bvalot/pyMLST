"""PyMLST entry commands and common parameters creation.

Subcommands are being instantiated dynamically from their respective folders.
"""

import os
import click

from click import Option

from pymlst import version
from pymlst.common import utils


class PyMlstCommand(click.MultiCommand):
    """Global PyMLST command."""

    def __init__(self, path, help_msg):
        """Initializes the command."""
        super().__init__(help='Subcommands are loaded from a '
                              'plugin folder dynamically')
        
        opt_help = dict(help_option_names=['-h', '--help'])

        opt_version = Option(['--version', '-v'], is_flag=True, callback=print_version,
                             expose_value=False, is_eager=True,
                             help='Prints PyMLST version.')
        self.params.append(opt_version)
        self.context_settings.update(opt_help)
        self.help = help_msg
        self.path = path

    def list_commands(self, ctx):
        """Lists the available commands.

        The commands are loaded dynamically from files within
        a directory.
        """
        cmd_names = []
        for filename in os.listdir(self.path):
            if filename.endswith('.py') and not filename.startswith('__init__'):
                cmd_names.append(filename[:-3])
        cmd_names.sort()
        return cmd_names

    def get_command(self, ctx, name):
        """Gets a command by name."""
        name_scope = {}
        cmd_file = os.path.join(self.path, name + '.py')
        try:
            with open(cmd_file) as file:
                code = compile(file.read(), cmd_file, 'exec')
                eval(code, name_scope, name_scope)
        except FileNotFoundError:
            raise click.ClickException(
                'Unknown sub-command \'{}\''.format(name))
        return name_scope['cli']


def print_version(ctx, param, value):
    """Prints the package version."""
    del param
    if not value or ctx.resilient_parsing:
        return
    click.echo('Version: ' + version.__version__)
    click.echo('Release: ' + version.__release__)
    ctx.exit()


py = PyMlstCommand(
    os.path.join(os.path.dirname(__file__), 'common', 'commands'),
    '''
    \b
    Common utility commands.
    
    Two pipelines are available:
      claMLST     for classical MLST analysis
      wgMLST      for Whole/Core genome MLST analysis''')

wg = PyMlstCommand(
    os.path.join(os.path.dirname(__file__), 'wg', 'commands'),
    'Whole/Core genome MLST commands.')

cla = PyMlstCommand(
    os.path.join(os.path.dirname(__file__), 'cla', 'commands'),
    'Classical MLST commands.')
