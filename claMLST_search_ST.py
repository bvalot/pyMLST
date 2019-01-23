#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Search ST number for an assembly"""

import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
import tempfile
import subprocess
import shutil

blat_exe = "blat"

desc = "Search ST number for a strain"
command = argparse.ArgumentParser(prog='claMLST_search_ST.py', \
    description=desc, usage='%(prog)s [options] genome database')
command.add_argument('-i', '--identity', nargs='?', \
    type=float, default=0.9, \
    help='Minimun identity to search gene (default=0.9)')
command.add_argument('-f', '--fasta', \
    type=argparse.FileType("w"), \
    help='Write fasta file with gene allele')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to BLAT executable (default=/usr/bin)')
command.add_argument('-o', '--output', default=sys.stdout, \
    type=argparse.FileType("w"), \
    help='Write ST search result to (default=stdout)')
command.add_argument('genome', \
    type=argparse.FileType("r"), \
    help='Genome of the strain')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database containing MLST shema')

def create_coregene(cursor, tmpfile):
    ref = int(1)
    cursor.execute('''SELECT DISTINCT gene FROM mlst''')
    all_rows = cursor.fetchall()
    coregenes = []
    for row in all_rows:
        cursor.execute('''SELECT sequence,gene FROM sequences WHERE allele=? and gene=?''', (1,row[0]))
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

def run_blat(path, genome, tmpfile, tmpout, identity):
    command = [path+blat_exe, '-maxIntron=4000', '-minIdentity='+str(identity*100),\
               genome.name, tmpfile.name, tmpout.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=sys.stderr)
    error = ""
    for line in iter(proc.stderr.readline,''):
        error += line
    if error != "":
        sys.stdout.write("Error during runnin BLAT\n")
        raise Exception(error)
    genes = {}
    for line in open(tmpout.name, 'r'):
        try:
            int(line.split()[0])
        except:
            continue
        psl = Psl(line)
        if psl.coverage == 1:
            genes.setdefault(psl.geneId(),[]).append(psl)
    if len(genes) == 0:
        raise Exception("No path found for the coregenome")
    return genes

def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs

class Psl:
    """A simple Psl class"""
    def __init__(self, pslline):
        pslelement = pslline.rstrip("\n").split("\t")
        if len(pslelement) != 21:
            raise Exception("Psl line have not 21 elements:\n"+pslline)
        self.pslelement = pslelement
        self.chro = pslelement[13]
        self.start = int(pslelement[15])
        self.end = int(pslelement[16])
        self.strand = pslelement[8]
        self.coverage = (float(self.pslelement[12]) - int(self.pslelement[11]))/int(self.pslelement[10])
        if self.coverage !=1 and self.coverage>=0.95:
            self.correct()
        
    def geneId(self):
        return self.pslelement[9]

    def correct(self):
        if int(self.pslelement[11]) != 0:
            diff = int(self.pslelement[11])
            if self.strand == "+":
                self.start = self.start - diff
            else:
                self.end = self.end + diff
        elif int(self.pslelement[10]) != int(self.pslelement[12]):
            diff = int(self.pslelement[10]) - int(self.pslelement[12])
            if self.strand == "+":
                self.end = self.end + diff
            else:
                self.start = self.start - diff   
        self.coverage = 1
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    genome = args.genome
    if args.identity<0 or args.identity > 1:
        raise Exception("Identity must be between 0 to 1")
    path = args.path.rstrip("/")+"/"
    if os.path.exists(path+blat_exe) is False:
        raise Exception("BLAT executable not found in folder: \n"+path)
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpout = tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=False)
    tmpout.close()
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()
        
        ##read coregene
        coregenes = create_coregene(cursor, tmpfile)
        tmpfile.close()

        ##BLAT analysis
        sys.stderr.write("Search coregene with BLAT\n")
        genes = run_blat(path, genome, tmpfile, tmpout, args.identity)
        sys.stderr.write("Finish run BLAT, found " + str(len(genes)) + " genes\n")
        
        ##Search sequence MLST
        seqs = read_genome(genome)
        sys.stderr.write("Search allele gene to database\n")
        # print(genes)
        allele = {i:[] for i in coregenes}
        st = {i:set() for i in coregenes}
        for coregene in coregenes:
            if coregene not in genes:
                allele.get(coregene).append("")
                continue
            for gene in genes.get(coregene):
                seq = seqs.get(gene.chro, None)
                if seq is None:
                    raise Exception("Chromosome ID not found " + gene.chro)
                ##get sequence
                sequence = str(seq.seq[gene.start:gene.end]).upper()
                if gene.strand !="+":
                    sequence = str(seq.seq[gene.start:gene.end].reverse_complement()).upper()
                ##verify complet sequence
                if len(sequence) != (gene.end-gene.start):
                    continue
                ##write fasta file with coregene
                if args.fasta is not None:
                    args.fasta.write(">"+coregene+"\n")
                    args.fasta.write(sequence+"\n")
                ##search allele
                cursor.execute('''SELECT allele FROM sequences WHERE sequence=? and gene=?''', \
                               (sequence, coregene))
                row = cursor.fetchone()
                if row is not None:
                    allele.get(coregene).append(str(row[0]))
                    cursor.execute('''SELECT st FROM mlst WHERE gene=? and allele=?''', \
                               (coregene,row[0]))
                    for row2 in cursor.fetchall():
                        st.get(coregene).add(row2[0])
                else:
                    allele.get(gene.geneId()).append("new")

        ##if only know allele or not found
        ##Seach st
        st_val = []
        if sum([len(i)==1 and i[0] != "new" for i in allele.values()]) == len(allele):
            tmp = None
            for s in st.values():
                if s:
                    if tmp is None:
                        tmp = s
                    else:
                        tmp = tmp.intersection(s)
            st_val = list(tmp)

        ##print result
        coregenes.sort()
        args.output.write("Sample\tST\t"+"\t".join(coregenes)+"\n")
        args.output.write(genome.name + "\t" + ";".join(map(str,st_val)))
        for coregene in coregenes:
            args.output.write("\t" + ";".join(map(str,allele.get(coregene))))
        args.output.write("\n")
        sys.stderr.write("FINISH\n")
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):        
            os.remove(tmpfile.name)
        if os.path.exists(tmpout.name):        
            os.remove(tmpout.name)
