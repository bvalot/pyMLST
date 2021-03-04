#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Remove strains to an wgMLST database"""
import os
import sys

import click
from pymlst.lib import sql
from pymlst.wg_commands.db.database import DatabaseCore
from pymlst.lib.benchmark import benchmark

desc = "Remove strain to a wgMLST database and sequences specificaly associated"


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of strains to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('strains',
                type=str, nargs=-1)
@benchmark
def cli(list, database, strains):
    """Remove strain to a wgMLST database and sequences specificaly associated"""

    if sql.ref in strains:
        raise Exception("Ref schema could not be remove from this database")

    # list strains to be removed

    all_strains = []
    if list is not None:
        for line in list.readlines():
            all_strains.append(line.rstrip("\n"))
    if strains is not None:
        for strain in strains:
            all_strains.append(strain)
    if len(all_strains) == 0:
        raise Exception("No strain to remove found.\n")
    all_strains = set(all_strains)

    try:
        db = DatabaseCore(os.path.abspath(database.name))

        for strain in all_strains:
            sys.stderr.write(strain + " : ")

            seqids = db.get_strain_sequences_ids(strain)
            if len(seqids) == 0:
                sys.stderr.write("Not found\n")
            else:
                sys.stderr.write("OK\n")

            db.remove_strain(strain)
            db.remove_orphan_sequences(seqids)

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
