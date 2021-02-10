#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import click
import os
from .__init__ import __version__

plugin_folder = os.path.join(os.path.dirname(__file__), 'wg_commands')


class CommandWG(click.MultiCommand):

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


cli = CommandWG(help='wgMLST\'s subcommands are loaded from a '
                     'plugin folder dynamically.')

if __name__ == '__main__':
    cli()
