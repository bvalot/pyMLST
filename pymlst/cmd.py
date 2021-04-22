#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import os

from click import Option

from pymlst import version
from pymlst.common import utils


class PyMlstCommand(click.MultiCommand):

    def __init__(self, path):
        super().__init__(help='Subcommands are loaded from a '
                              'plugin folder dynamically')
        opt_version = Option(['--version', '-v'], is_flag=True, callback=print_version,
                             expose_value=False, is_eager=True,
                             help='Prints PyMLST version.')
        opt_debug = Option(['--debug', '-d'], is_flag=True, callback=set_debug,
                             expose_value=False, is_eager=False,
                             help='Sets the debug mode ON.')
        self.params.append(opt_version)
        self.params.append(opt_debug)
        self.path = path

    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir(self.path):
            if filename.endswith('.py') and not filename.startswith('__init__'):
                rv.append(filename[:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        ns = {}
        fn = os.path.join(self.path, name + '.py')
        with open(fn) as f:
            code = compile(f.read(), fn, 'exec')
            eval(code, ns, ns)
        return ns['cli']


def set_debug(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    utils.create_logger(verbose=True)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('Version: ' + version.__version__)
    click.echo('Release: ' + version.__release__)
    ctx.exit()


py = PyMlstCommand(os.path.join(os.path.dirname(__file__), 'common', 'commands'))
wg = PyMlstCommand(os.path.join(os.path.dirname(__file__), 'wg', 'commands'))
cla = PyMlstCommand(os.path.join(os.path.dirname(__file__), 'cla', 'commands'))
