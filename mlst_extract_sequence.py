#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Get sequence from wgMLST database"""

import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
import tempfile
import subprocess
from lib import __version__
import lib.sql as sql

mafft_exe = "mafft"

desc = "Get sequences from wgMLST database"
command = argparse.ArgumentParser(prog='mlst_extract_sequence.py', \
    description=desc, usage='%(prog)s [options] database')
command.add_argument('-o', '--output', nargs='?', \
    type=argparse.FileType("w"), default=sys.stdout, \
    help='Output result on fasta format in (default:stdout)')
command.add_argument('-l', '--liste', nargs='?', \
    type=argparse.FileType("r"), default=None, \
    help='List of coregene to extract (default:all)')
command.add_argument('-a', '--align', action='store_true', \
    help='Report a concatened multi-fasta file instead of only gene files (default:No)')
command.add_argument('-r', '--realign', action='store_true', \
    help='Realign gene with same length (Default:No)')
command.add_argument('-m', '--mincover', nargs='?', \
    type=int, default=1, \
    help='Minimun number of strain found to keep a coregene (default:1)')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to mafft executable (default=/usr/bin)')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to store MLST')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)

def run_mafft(path, tmpfile):
    command = [path+mafft_exe, '--quiet', tmpfile.name]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    genes = {}
    ids = None
    seq=""
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
            raise Exception("Problem during run mafft" + line)
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
    """Add sequence to multialign, take first gene in case of repeat"""
    size = 0
    if len(seqs)>0:
        size = len(seqs[0][2])
    for s in strains:
        seq = [i[2] for i in seqs if s in i[1]]
        if len(seq) == 0:
            sequences.get(s).append('-' * size)
        elif len(seq) == 1:
            sequences.get(s).append(seq[0])
        else:
            raise Exception("repeat gene must be excluded for align export\n")


if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    output = args.output
    path = args.path.rstrip("/")+"/"
    if os.path.exists(path+mafft_exe) is False:
        raise Exception("Mafft executable not found in folder: \n"+path)

    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpfile.close()
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        ##index old database
        sql.index_database(cursor)
        
        ##Minimun number of strain
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (sql.ref,))
        strains = [i[0] for i in cursor.fetchall()]
        if args.mincover < 1 or args.mincover > len(strains):
            raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))

        ## Coregene
        coregene = []
        if args.liste is not None:
            coregene = [l.rstrip("\n") for l in iter(args.liste.readline, '')]
        else:
            cursor.execute('''SELECT gene, COUNT(distinct souche) FROM mlst
                              WHERE souche != ?
                              GROUP BY gene''', (sql.ref,))
            coregene = [l[0] for l in cursor.fetchall() if l[1] >= args.mincover]

        if args.align is False:
            ##no multialign
            for g in coregene:
                seqs = get_sequences_for_gene(cursor, g)
                for seq in seqs:
                    output.write(">"+ g + "|" + str(seq[0]) + " " \
                                 + ";".join(seq[1]) + "\n")
                    output.write(seq[2] + "\n")
        else:
            ## multialign
            ## search duplicate
            dupli = sql.get_duplicate_gene(cursor)
            
            sequences = {s:[] for s in strains}
            for i,g in enumerate(coregene):
                sys.stderr.write(str(i+1) + "/" + str(len(coregene)) + " | " + g + "     ")
                ##geneid = get_geneids(cursor, g)
                if g in dupli:
                    sys.stderr.write("No: Repeat gene\n")
                    continue
                seqs = get_sequences_for_gene(cursor, g)
                size = set()
                for seq in seqs:
                    size.add(len(seq[2]))
                if len(size) == 1 and args.realign is False:
                    sys.stderr.write("Direct")
                    add_sequence_strain(seqs, strains, sequences)
                else:
                    sys.stderr.write("Align")
                    write_tmp_seqs(tmpfile, seqs)
                    corrseqs = run_mafft(path, tmpfile)
                    for seq in seqs:
                        seq[2] = corrseqs.get(seq[0])
                    add_sequence_strain(seqs, strains, sequences)
                sys.stderr.write("\n")

            ##output align result
            for s in strains:
                output.write('>'+ s + "\n")
                output.write("\n".join(sequences.get(s)) + "\n")
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):        
            os.remove(tmpfile.name)
