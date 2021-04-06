#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL


import sys
import click

from pymlst.wg.core import find_subgraph


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output group files (default:stdout)')
@click.option('--threshold', '-t',
              type=int, default=50,
              help='Minimum distance to conserve '
                   'for extraction of group (default:50)')
@click.option('--count', '-c',
              type=click.File('w'),
              help='Output count file')
@click.argument('distance',
                type=click.File('r'))
def cli(output, threshold, count, distance):
    """Search group os strain at a distance treshold"""

    find_subgraph(threshold, count, distance, output)
