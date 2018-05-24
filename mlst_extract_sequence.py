#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Get sequence from wgMLST database"""

import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
import tempfile
import subprocess

mafft_exe = "mafft"
ref = 'ref'

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
    help='Minimun number of strain found to conserved a coregene (default:1)')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to mafft executable (default=/usr/bin)')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to stock MLST')


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
            genes[int(ids)] = seq
            ids = line.lstrip('>').rstrip("\n")
            seq = ""
        elif ids is not None:
            seq += line.rstrip("\n")
        else:
            raise Exception("Problem during run mafft" + line)
    if seq != "":
        genes[int(ids)] = seq     
    # print(genes)
    return genes

# cursor.execute("SELECT id,sequence FROM sequences WHERE id in (%s)" % ",".join(map(str,geneid)))

def get_geneids(cursor, gene):
    cursor.execute('''SELECT souche, seqid FROM mlst WHERE gene=? and souche!=?''', (gene,ref))
    geneids = {}
    for row in cursor.fetchall():
        geneids.setdefault(row[0], []).append(row[1])
    return geneids

def get_sequence(cursor, ids):
    cursor2.execute("SELECT id,sequence FROM sequences WHERE id in (%s)" % ",".join(map(str,ids)))
    seqs = {}
    for row in cursor.fetchall():
        seqs[row[0]] = row[1]
    return seqs

def write_tmp_seqs(tmpfile, seqs):
    tmp = open(tmpfile.name, 'w+t')
    for s in seqs.items():
        tmp.write(">"+str(s[0])+"\n"+s[1]+"\n")
    tmp.close()

def add_sequence_strain(geneid, seqs, strains, sequences):
    """Add sequence to multialign, take first gene in case of repeat""" 
    for s in strains:
        if not geneid.get(s):
            sequences.get(s).append('-' * len(seqs.values()[0]))
        else:
            sequences.get(s).append(seqs.get(geneid.get(s)[0]))


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

        coregene = []
        if args.liste is not None:
            coregene = [l.rstrip("\n") for l in iter(args.liste.readline, '')]
        else:
            cursor.execute('''SELECT DISTINCT gene FROM mlst''')
            coregene = [l[0] for l in cursor.fetchall()]

        ##Minimun number of strain
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (ref,))
        strains = [i[0] for i in cursor.fetchall()]
        if args.mincover <= 1 or args.mincover > len(strains):
            raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))
        
        if args.align is False:
            ##no multialign
            for g in coregene:
                geneid = get_geneids(cursor, g)
                if len(geneid) < args.mincover:
                    continue
                seqs = get_sequence(cursor2, set([item for sublist in geneid.values() for item in sublist]))
                for seq in seqs.items():
                    output.write(">"+ g + "|" + str(seq[0]) + " ")
                    output.write(";".join([i[0] for i in geneid.items() if seq[0] in i[1]]) + "\n")
                    output.write(seq[1] + "\n")
        else:
            ##multialign
            sequences = {s:[] for s in strains}
            for i,g in enumerate(coregene):
                sys.stderr.write(str(i+1) + "/" + str(len(coregene)) + " | " + g + "     ")
                geneid = get_geneids(cursor, g)
                if len(geneid) < args.mincover :
                    sys.stderr.write("No: Less mincover \n")
                    continue
                if max(map(len,geneid.values()))>1:
                    sys.stderr.write("No: Repeat gene\n")
                    continue
                seqs = get_sequence(cursor2, set([item for sublist in geneid.values() for item in sublist]))
                if len(set(map(len,seqs.values()))) == 1 and args.realign is False:
                    sys.stderr.write("Direct")
                    add_sequence_strain(geneid, seqs, strains, sequences)
                else:
                    sys.stderr.write("Align")
                    write_tmp_seqs(tmpfile, seqs)
                    corrseqs = run_mafft(path, tmpfile)
                    add_sequence_strain(geneid, corrseqs, strains, sequences)
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
