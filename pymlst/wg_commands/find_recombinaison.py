#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Search potential recombinaison from wgMLST database export"""

import sys
import click

from pymlst.api.wgmlst import find_recombination

desc = "Search potential gene recombinaison from wgMLST database export"


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output number of variations by genes (default:stdout)')
@click.argument('genes',
                type=click.File('r'))
@click.argument('alignment',
                type=click.File('r'))
def cli(output, genes, alignment):
    """Search potential gene recombinaison from wgMLST database export"""

    find_recombination(genes, alignment, output)