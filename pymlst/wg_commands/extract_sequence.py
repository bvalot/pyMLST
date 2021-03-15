#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Get a sequence from a wgMLST database"""

import sys
import os

import click

from pymlst.api.core import open_wg
from pymlst.lib.benchmark import benchmark
from pymlst.api.extractors import SequenceExtractor

mafft_exe = '/usr/bin/mafft'

desc = "Get a sequence from a wgMLST database"

@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output result in fasta format (default:stdout)')
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='List of coregenes to extract (default:all)')
@click.option('--align', '-a',
              is_flag=True,
              help='Report a concatened multi-fasta file '
              'instead of only gene files (default:No)')
@click.option('--realign', '-r',
              is_flag=True,
              help='Realign genes with same length (Default:No)')
@click.option('--mincover', '-m',
              type=int, default=1,
              help='Minimun number of strain found '
              'to keep a coregene (default:1)')
@click.argument('database',
                type=click.File('r'))
@benchmark  # TEST
def cli(output, list, align, realign, mincover, database):
    """Get sequences from a wgMLST database"""

    database.close()

    with open_wg(os.path.abspath(database.name)) as mlst:
        mlst.extract(SequenceExtractor(list, align, realign, mincover), output)
