#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Extract MLST table from an wgMLST database"""
import os
import sys

import click

from pymlst.api.extractors import TableExtractor, ExportType
from pymlst.api.core import open_wg

from pymlst.lib.benchmark import benchmark

desc = "Extract MLST table from a wgMLST database"


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Export MLST table to (default=stdout)')
@click.option('--export', '-e',
              type=click.Choice(ExportType.list_types()),
              default='mlst',
              help='Defined the export format')
@click.option('--count', '-c',
              is_flag=True,
              help='In strain mode, count the number of gene present in the database')
@click.option('--mincover', '-m',
              type=int, default=0,
              help='Minimun number of strain found to keep a gene (default:0)')
@click.option('--keep', '-k',
              is_flag=True,
              help='Keep only gene with different allele (omit missing)')
@click.option('--duplicate', '-d',
              is_flag=True, default=True,
              help='Conserve duplicate gene (default remove)')
@click.option('--inverse', '-V',
              is_flag=True,
              help='Keep only gene that do not ' \
              'meet the filter of mincover or keep options')
@click.argument('database', type=click.File('r'))
@benchmark
def cli(output, export, count, mincover, keep, duplicate, inverse, database):
    """Extract MLST table from a wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.extract(TableExtractor(export, count, mincover, keep, duplicate, inverse), output)
