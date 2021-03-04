#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Remove gene to a wgMLST database"""
import os
import sys

import click

from pymlst.wg_commands.db.database import DatabaseCore

desc = "Remove gene to a wgMLST database and sequences specificaly associated"


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of genes to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('genes',
                type=str, nargs=-1)
def cli(list, genes, database):
    """Remove gene to a wgMLST database"""

    # list genes to removed

    all_genes = []
    if list is not None:
        for line in list.readlines():
            all_genes.append(line.rstrip("\n"))
    if genes is not None:
        for gene in genes:
            all_genes.append(gene)
    if len(all_genes) == 0:
        raise Exception("No gene to remove found.\n")
    all_genes = set(all_genes)

    try:
        db = DatabaseCore(os.path.abspath(database.name))

        for gene in all_genes:
            sys.stderr.write(gene + " : ")

            seqids = db.get_gene_sequences_ids(gene)
            if len(seqids) == 0:
                sys.stderr.write("Not found\n")
            else:
                sys.stderr.write("OK\n")

            db.remove_gene(gene)
            db.remove_orphan_sequences(seqids)

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
