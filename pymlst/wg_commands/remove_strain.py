#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Remove strains to an wgMLST database"""

import sys
import sqlite3
import click
from pymlst.lib import sql

desc = "Remove strain to a wgMLST database and sequences specificaly associated"


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of strains to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('strains',
                type=str, nargs=-1)
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
        db = sqlite3.connect(database.name)
        cursor = db.cursor()

        # index old database
        sql.index_database(cursor)

        for strain in all_strains:
            sys.stderr.write(strain + "     ")

            # Search seq ids
            cursor.execute('''SELECT seqid FROM mlst WHERE souche=?''', (strain,))
            seqids = cursor.fetchall()
            if len(seqids) == 0:
                raise Exception("Strain name not found in database\n" + strain)

            # remove sample
            cursor.execute('''DELETE FROM mlst WHERE souche=? ''', (strain,))

            # remove seqs if no other strain have this seq
            for seqid in seqids:
                cursor.execute('''DELETE from sequences as s
                                  where not exists (
                                  select 1 from mlst where seqid=s.id)
                                  and id=? ''', (seqid[0],))
            sys.stderr.write("OK\n")

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
