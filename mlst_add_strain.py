#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Add a strain to wgMLST database"""

import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
import tempfile
import subprocess
import shutil

gmapbuild_exe = "gmap_build"
gmap_exe = "gmap"

desc = "Add a strain to the wgMLST database"
command = argparse.ArgumentParser(prog='mlst_add_strain.py', \
    description=desc, usage='%(prog)s [options] genome database')
command.add_argument('-s', '--strain', nargs='?', \
    type=str, default=None, \
    help='Name of the strain (default:genome name)')
command.add_argument('-i', '--identity', nargs='?', \
    type=float, default=0.95, \
    help='Minimun identity to search gene (default=0.95)')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to gmap executable (default=/usr/bin)')
command.add_argument('genome', \
    type=argparse.FileType("r"), \
    help='Genome of the strain')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to stock MLST')

def create_coregene(cursor, tmpfile):
    ref = "ref"
    cursor.execute('''SELECT gene, seqid FROM mlst WHERE souche=?''', (ref,))
    all_rows = cursor.fetchall()
    coregenes = []
    for row in all_rows:
        cursor.execute('''SELECT sequence FROM sequences WHERE id=?''', (row[1],))
        tmpfile.write('>' + row[0] + "\n" + cursor.fetchone()[0] + "\n")
        coregenes.append(row[0])
    return coregenes

def insert_sequence(cursor, sequence):
    try:
        cursor.execute('''INSERT INTO sequences(sequence) VALUES(?)''', (sequence,))
        return cursor.lastrowid
    except sqlite3.IntegrityError:
        cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (sequence,))
        return cursor.fetchone()[0]

def build_gmap(path, genome, name):
    command = [path+gmapbuild_exe, '-k', '15', '-D', tempfile.gettempdir(), '-d', name, genome.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=sys.stdout)
    error = ""
    for line in iter(proc.stderr.readline,''):
        error += line
    # if error != "":
    #     sys.stdout.write("Error during build gmap\n")
    #     raise Exception(error)

def run_gmap(path, name, tmpfile, identity):
    command = [path+gmap_exe, '--gff3-add-separators', '0', '--min-identity', str(identity),\
               '--min-trimmed-coverage', '1', '-f', '3', '--nosplicing',\
               '-D', tempfile.gettempdir(), '-d', name, tmpfile.name]
    proc = subprocess.Popen(command, stderr=sys.stdout, stdout=subprocess.PIPE)
    genes = {}
    for line in iter(proc.stdout.readline, ''):
        if line[0] =="#":
            continue
        gff = Gff(line)
        genes.setdefault(gff.geneId(),[]).append(gff)
    if len(genes) == 0:
        raise Exception("No path found for the coregenome")
    return genes

def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs


class Gff:
    """A simple Gff class"""
    def __init__(self, gffline):
        gffelement = gffline.strip().split("\t")
        if len(gffelement) != 9:
            raise Exception("Gff line have not 9 elements:\n"+gffline)
        self.chro = gffelement[0]
        self.start = int(gffelement[3])
        self.end = int(gffelement[4])
        self.identity = int(gffelement[5])
        self.values = {a.split("=")[0]:a.split("=")[1] for a in gffelement[8].split(";")\
                       if len(a.split("="))==2}
        self.strand = gffelement[6]
        
    def is_in(self, chro, pos):
        if chro != self.chro:
            return False
        elif pos>self.start and pos<self.end:
            return True
        else:
            return False

    def geneId(self):
        return self.values.get("Name", None)
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    genome = args.genome
    name = args.strain
    if args.identity<0 or args.identity > 1:
        raise Exception("Identity must be between 0 to 1")
    path = args.path.rstrip("/")+"/"
    if os.path.exists(path+gmap_exe) is False:
        raise Exception("Gmap executable not found in folder: \n"+path)
    if name is None:
        name = genome.name.split('/')[-1]
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)

    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        ##verify strain not alrealdy on the database
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche=?''', (name,))
        if cursor.fetchone() is not None:
            raise Exception("Strain is already present in database:\n"+name)
        
        ##read coregene
        coregenes = create_coregene(cursor, tmpfile)
        tmpfile.close()

        ##gmap analysis
        sys.stderr.write("Build GMAP database      ")
        build_gmap(path, genome, name)
        sys.stderr.write("OK\n")

        sys.stderr.write("Search coregene with GMAP      ")
        genes = run_gmap(path, name, tmpfile, args.identity)
        sys.stderr.write("OK\n")
        sys.stderr.write("Finish run GMAP, found " + str(len(genes)) + " genes\n")
        
        ##add sequence MLST
        seqs = read_genome(genome)
        sys.stderr.write("Add new MLST gene to database\n")
        for coregene in coregenes:
            if coregene not in genes:
                continue
            gene = genes.get(coregene)
            if len(gene) > 1 :
                sys.stderr.write("Gene found in duplicate: " + coregene + " \n")
                continue
            seq = seqs.get(gene[0].chro, None)
            if seq is None:
                raise Exception("Chromosome ID not found " + gene[0].chro)
            ##add sequence and MLST
            if gene[0].strand =="+":
                seqid = insert_sequence(cursor, str(seq.seq[gene[0].start-1:gene[0].end]))
            else:
                seqid = insert_sequence(cursor, str(seq.seq[gene[0].start-1:gene[0].end].reverse_complement()))
            cursor2.execute('''INSERT INTO mlst(souche, gene, seqid)
                              VALUES(?,?,?)''', (name, gene[0].geneId(), seqid))
        db.commit()
        sys.stderr.write("FINISH\n")
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()
        if os.path.exists(tmpfile.name):        
            os.remove(tmpfile.name)
        gmapdir = tempfile.gettempdir()+"/"+name
        if os.path.exists(gmapdir):
            shutil.rmtree(gmapdir)
