#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import os

from click import Option

from pymlst import version

plugin_folder = os.path.join(os.path.dirname(__file__), 'cla_commands')


class CommandCLA(click.MultiCommand):

    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir(plugin_folder):
            if filename.endswith('.py'):
                rv.append(filename[:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        ns = {}
        fn = os.path.join(plugin_folder, name + '.py')
        with open(fn) as f:
            code = compile(f.read(), fn, 'exec')
            eval(code, ns, ns)
        return ns['cli']


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('Version: ' + version.__version__)
    click.echo('Release: ' + version.__release__)
    ctx.exit()


cli = CommandCLA(help='claMLST\'s subcommands are loaded from a '
                      'plugin folder dynamically.')

opt = Option(['--version', '-v'], is_flag=True, callback=print_version,
             expose_value=False, is_eager=True,
             help='Prints PyMLST version.')

cli.params.append(opt)

if __name__ == '__main__':
    cli()
