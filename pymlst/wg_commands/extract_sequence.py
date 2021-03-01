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

import tempfile
import subprocess

from pymlst.lib import sql
from pymlst.lib.benchmark import benchmark
from pymlst.wg_commands.db.database import DatabaseCore

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
            raise Exception("A problem occurred while running mafft" + str(line))
    if seq != "":
        genes[int(ids)] = seq.upper()

    return genes


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
@benchmark  # TEST
def cli(output, list, align, realign, mincover, database):
    """Get sequences from a wgMLST database"""

    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpfile.close()

    try:
        db = DatabaseCore(os.path.abspath(database.name))

        # Minimun number of strain
        strains = [i[0] for i in db.get_different_souches(sql.ref)]
        if mincover < 1 or mincover > len(strains):
            raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))

        #  Coregene
        coregene = []
        if list is not None:
            coregene = [l.rstrip("\n") for l in iter(list.readline, '')]
        else:
            coregene = [l[0] for l in db.get_genes_coverages(sql.ref) if l[1] >= mincover]

        sys.stderr.write("Number of gene to analyse : " + str(len(coregene)) + "\n")

        if align is False:
            # no multialign
            for g in coregene:
                seqs = db.get_gene_sequences(g, sql.ref)
                for seq in seqs:
                    output.write(">" + g + "|" + str(seq[0]) + " "
                                 + ";".join(seq[1]) + "\n")
                    output.write(seq[2] + "\n")
        else:
            # multialign
            # search duplicate
            dupli = db.get_duplicated_genes(sql.ref)
            for d in dupli:
                print('Duplicated: ', d)

            sequences = {s:[] for s in strains}
            for i, g in enumerate(coregene):
                sys.stderr.write(str(i+1) + "/" + str(len(coregene)) + " | " + g + "     ")
                if g in dupli:
                    sys.stderr.write("No: Repeat gene\n")
                    continue
                seqs = db.get_gene_sequences(g, sql.ref)
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

            db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):
            os.remove(tmpfile.name)
