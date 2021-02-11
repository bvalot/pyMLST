#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Get a sequence from a wgMLST database"""

import sys
import os

import sqlite3
import click

import tempfile
import subprocess

from pymlst.lib import sql

mafft_exe = '/usr/bin/mafft'

desc = "Get a sequence from a wgMLST database"


def run_mafft(path, tmpfile):
    command = [path, '--quiet', tmpfile.name]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    genes = {}
    ids = None
    seq = ""
    for line in iter(proc.stdout.readline, ''):
        if ids is None and line[0] == '>':
            ids = line.lstrip('>').rstrip("\n")
        elif ids is not None and line[0] == '>':
            genes[int(ids)] = seq.upper()
            ids = line.lstrip('>').rstrip("\n")
            seq = ""
        elif ids is not None:
            seq += line.rstrip("\n")
        else:
            raise Exception("A problem occurred while running mafft" + line)
    if seq != "":
        genes[int(ids)] = seq.upper()
    # print(genes)
    return genes


def get_sequences_for_gene(cursor, gene):
    cursor.execute('''SELECT m.seqid, group_concat(m.souche, ";"), s.sequence
                      FROM mlst as m
                      JOIN sequences as s
                      ON m.seqid = s.id
                      WHERE m.souche != ?
                      AND m.gene = ?
                      GROUP BY m.seqid''', (sql.ref, gene))
    seqs = []
    for seq in cursor.fetchall():
        tmp = seq[1].split(";")
        tmp.sort()
        seqs.append([seq[0], tmp, seq[2]])
    return seqs


def write_tmp_seqs(tmpfile, seqs):
    tmp = open(tmpfile.name, 'w+t')
    for s in seqs:
        tmp.write(">"+str(s[0])+"\n"+s[2]+"\n")
    tmp.close()


def add_sequence_strain(seqs, strains, sequences):
    """Add a sequence to multi-align, take the first gene in case of repetition"""
    size = 0
    if len(seqs) > 0:
        size = len(seqs[0][2])
    for s in strains:
        seq = [i[2] for i in seqs if s in i[1]]
        if len(seq) == 0:
            sequences.get(s).append('-' * size)
        elif len(seq) == 1:
            sequences.get(s).append(seq[0])
        else:
            raise Exception("repeated gene must be excluded in order to align export\n")


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
def cli(output, list, align, realign, mincover, database):
    """Get sequences from a wgMLST database"""

    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpfile.close()

    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        # index old database
        sql.index_database(cursor)

        # Minimun number of strain
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (sql.ref,))
        strains = [i[0] for i in cursor.fetchall()]
        if mincover < 1 or mincover > len(strains):
            raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))

        #  Coregene
        coregene = []
        if list is not None:
            coregene = [l.rstrip("\n") for l in iter(list.readline, '')]
        else:
            cursor.execute('''SELECT gene, COUNT(distinct souche) FROM mlst
                              WHERE souche != ?
                              GROUP BY gene''', (sql.ref,))
            coregene = [l[0] for l in cursor.fetchall() if l[1] >= mincover]

        sys.stderr.write("Number of gene to analyse : " + str(len(coregene)) + "\n")

        if align is False:
            # no multialign
            for g in coregene:
                seqs = get_sequences_for_gene(cursor, g)
                for seq in seqs:
                    output.write(">" + g + "|" + str(seq[0]) + " "
                                 + ";".join(seq[1]) + "\n")
                    output.write(seq[2] + "\n")
        else:
            # multialign
            # search duplicate
            dupli = sql.get_duplicate_gene(cursor)

            sequences = {s:[] for s in strains}
            for i, g in enumerate(coregene):
                sys.stderr.write(str(i+1) + "/" + str(len(coregene)) + " | " + g + "     ")
                # geneid = get_geneids(cursor, g)
                if g in dupli:
                    sys.stderr.write("No: Repeat gene\n")
                    continue
                seqs = get_sequences_for_gene(cursor, g)
                size = set()
                for seq in seqs:
                    size.add(len(seq[2]))
                if len(size) == 1 and realign is False:
                    sys.stderr.write("Direct")
                    add_sequence_strain(seqs, strains, sequences)
                else:
                    sys.stderr.write("Align")
                    write_tmp_seqs(tmpfile, seqs)
                    corrseqs = run_mafft(mafft_exe, tmpfile)
                    for seq in seqs:
                        seq[2] = corrseqs.get(seq[0])
                    add_sequence_strain(seqs, strains, sequences)
                sys.stderr.write("\n")

            # output align result
            for s in strains:
                output.write('>' + s + "\n")
                output.write("\n".join(sequences.get(s)) + "\n")
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):
            os.remove(tmpfile.name)
