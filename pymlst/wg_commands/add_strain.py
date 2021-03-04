#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Add a strain to the wgMLST database"""

import sys
import os

from Bio import SeqIO
import click
import time
from pymlst.lib import psl
from pymlst.lib import blat
from pymlst.lib import sql
from pymlst.wg_commands.db.database import DatabaseWG

blat_path = '/usr/bin/'

desc = "Add a strain to the wgMLST database"


def create_coregene(db, tmpfile):
    ref_genes = db.get_gene_by_souche(sql.ref)
    coregenes = []
    for row in ref_genes:
        tmpfile.write('>' + row.gene + "\n" + row.sequence + "\n")
        coregenes.append(row[0])
    return coregenes


def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs


@click.command()
@click.option('--strain', '-s',
              type=str, default=None,
              help='Name of the strain (default:genome name)')
@click.option('--identity', '-i',
              type=float, default=0.95,
              help='Minimum identity to search gene (default=0.95)')
@click.option('--coverage', '-c',
              type=float, default=0.9,
              help='Minimum coverage to search gene (default=0.9)')
@click.argument('genome',
                type=click.File("r"))
@click.argument('database',
                type=click.File("r"))
def cli(strain, identity, coverage, genome, database):
    """Add a strain to the wgMLST database"""
    start = time.time()

    if identity < 0 or identity > 1:
        raise Exception("Identity must be between 0 and 1")
    path = blat.test_blat_exe(blat_path)

    name = strain
    if name is None:
        name = genome.name.split('/')[-1]
    if ";" in name:
        raise Exception("Strain name must not contains special ';'\n")

    tmpfile, tmpout = blat.blat_tmp()

    try:
        db = DatabaseWG(os.path.abspath(database.name))

        # verify that the strain is not already in the database
        if db.contains_souche(name):
            raise Exception("Strain is already present in database:\n"+name)

        # read coregene
        coregenes = create_coregene(db, tmpfile)
        tmpfile.close()

        # BLAT analysis
        sys.stderr.write("Search coregene with BLAT\n")
        genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage)
        sys.stderr.write("Finish run BLAT, found " + str(len(genes)) + " genes\n")

        # add MLST sequence
        seqs = read_genome(genome)
        bad = 0
        for coregene in coregenes:
            if coregene not in genes:
                continue
            for gene in genes.get(coregene):
                seq = seqs.get(gene.chro, None)
                if seq is None:
                    raise Exception("Chromosome ID not found " + gene.chro)

                # Correct coverage
                if gene.coverage != 1:
                    if gene.searchPartialCDS(seq, coverage) is False:
                        sys.stderr.write("Gene " + gene.geneId() + " partial: removed\n")
                        bad += 1
                        continue
                    else:
                        sys.stderr.write("Gene " + gene.geneId() + " fill: added\n")

                # Verify CDS
                if psl.testCDS(gene.getSequence(seq), False) is False:
                    if gene.searchCorrectCDS(seq, coverage) is False:
                        sys.stderr.write("Gene " + gene.geneId() + " not correct: removed\n")
                        bad += 1
                        continue
                    else:
                        sys.stderr.write("Gene " + gene.geneId() + " correct: added\n")

                # add sequence and MLST
                sequence = gene.getSequence(seq)

                # Insert data in database
                seqid = db.add_sequence(str(sequence))[1]
                db.add_mlst(name, gene.geneId(), seqid)

        db.commit()
        sys.stderr.write("Add " + str(len(genes)-bad) + " new MLST gene to database\n")
        sys.stderr.write("FINISH\n")

        end = time.time()
        print('Elapsed time: ', end - start)
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):
            os.remove(tmpfile.name)
        if os.path.exists(tmpout.name):
            os.remove(tmpout.name)
