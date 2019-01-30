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

blat_exe = "blat"

desc = "Add a strain to the wgMLST database"
command = argparse.ArgumentParser(prog='mlst_add_strain.py', \
    description=desc, usage='%(prog)s [options] genome database')
command.add_argument('-s', '--strain', nargs='?', \
    type=str, default=None, \
    help='Name of the strain (default:genome name)')
command.add_argument('-i', '--identity', nargs='?', \
    type=float, default=0.95, \
    help='Minimun identity to search gene (default=0.95)')
command.add_argument('-c', '--coverage', nargs='?', \
    type=float, default=0.9, \
    help='Minimun coverage to search gene (default=0.9)')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to BLAT executable (default=/usr/bin)')
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
        seq = cursor.fetchone()[0]
        tmpfile.write('>' + row[0] + "\n" + seq + "\n")
        coregenes.append((row[0], seq))
    return coregenes

def insert_sequence(cursor, sequence):
    try:
        cursor.execute('''INSERT INTO sequences(sequence) VALUES(?)''', (sequence,))
        return cursor.lastrowid
    except sqlite3.IntegrityError:
        cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (sequence,))
        return cursor.fetchone()[0]

def run_blat(path, name, genome, tmpfile, tmpout, identity, coverage):
    command = [path+blat_exe, '-maxIntron=20', '-fine', '-minIdentity='+str(identity*100),\
               genome.name, tmpfile.name, tmpout.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=sys.stdout)
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
        if psl.coverage >=coverage and psl.coverage <= 1:
            genes.setdefault(psl.geneId(),[]).append(psl)
    if len(genes) == 0:
        raise Exception("No path found for the coregenome")
    return genes

def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs

def testCDS(seq, reverse):
    try:
        if reverse:
            seq.reverse_complement().translate(table="Bacterial", cds=True)
        else:
            seq.translate(table="Bacterial", cds=True)
    except:
        return False
    return True

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
        self.rstart = int(pslelement[11])
        self.rend = int(pslelement[12])
        self.rtotal = int(pslelement[10])
        self.coverage = (float(self.rend) - self.rstart)/self.rtotal
        # if self.coverage !=1 and self.coverage>=0.95:
        #     self.correct()
        
    def geneId(self):
        return self.pslelement[9]

    def searchCDS(self, seq):
        ##modifs start and stop not create
        if self.rstart !=0 and self.rend != self.rtotal:
            return False
        ##modifs start
        if self.rstart !=0:
            modulo = (self.end-self.start)%3
            diff = self.rstart
            if self.strand == "+":
                val = [i for i in range(self.start+modulo, self.start-diff*2, -3) \
                       if testCDS(seq.seq[i:self.end], False)]
                if len(val)==1:
                    self.start=val[0]
                    return True
                elif len(val) >1:
                    sys.stderr.write("Choice best start for gene " + gene.geneId() + "\n" + str(val) + "\n")
                    self.start = val[0]
                    return True
                else:
                    return False
            else:
                val = [i for i in range(self.end-modulo, self.end+diff*2, 3) \
                       if testCDS(seq.seq[self.start:i], True)]
                if len(val) == 1:
                    self.end = val[0]
                    return True
                elif len(val) >1:
                    sys.stderr.write("Choice best start for gene " + gene.geneId() + "\n" + str(val) + "\n")
                    self.end = val[0]
                    return True
                else:
                    return False
        ##modifs end
        elif self.rend != self.rtotal:
            modulo = (self.end-self.start)%3
            diff = self.rtotal - self.rend
            if self.strand == "+":
                val = [i for i in range(self.end-modulo, self.end+diff*2, 3) \
                       if testCDS(seq.seq[self.start:i], False)]
                if len(val) == 1:
                    self.end = val[0]
                    return True
                else:
                    return False
            else:
                val = [i for i in range(self.start+modulo, self.start-diff*2, -3) \
                       if testCDS(seq.seq[i:self.end], True)]
                if len(val) == 1:
                    self.start = val[0]
                    return True
                else:
                    return False
        else:
            raise Exception("A problem of coverage for gene " + gene.geneId())
    
    # def correct(self):
    #     if int(self.pslelement[11]) != 0:
    #         diff = int(self.pslelement[11])
    #         if self.strand == "+":
    #             self.start = self.start - diff
    #         else:
    #             self.end = self.end + diff
    #     elif int(self.pslelement[10]) != int(self.pslelement[12]):
    #         diff = int(self.pslelement[10]) - int(self.pslelement[12])
    #         if self.strand == "+":
    #             self.end = self.end + diff
    #         else:
    #             self.start = self.start - diff   
    #     self.coverage = 1
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    genome = args.genome
    name = args.strain
    if args.identity<0 or args.identity > 1:
        raise Exception("Identity must be between 0 to 1")
    path = args.path
    if path:
        path = path.rstrip("/")+"/"
        if os.path.exists(path+blat_exe) is False:
            raise Exception("BLAT executable not found in folder: \n"+path)
    if name is None:
        name = genome.name.split('/')[-1]
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpout = tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=False)
    tmpout.close()
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        ##verify strain not already on the database
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche=?''', (name,))
        if cursor.fetchone() is not None:
            raise Exception("Strain is already present in database:\n"+name)
        
        ##read coregene
        coregenes = create_coregene(cursor, tmpfile)
        tmpfile.close()

        ##BLAT analysis
        sys.stderr.write("Search coregene with BLAT\n")
        genes = run_blat(path, name, genome, tmpfile, tmpout, args.identity, args.coverage)
        sys.stderr.write("Finish run BLAT, found " + str(len(genes)) + " genes\n")
        
        ##add sequence MLST
        seqs = read_genome(genome)
        bad = 0
        for coregene in coregenes:
            if coregene[0] not in genes:
                continue
            for gene in genes.get(coregene[0]):
                seq = seqs.get(gene.chro, None)
                if seq is None:
                    raise Exception("Chromosome ID not found " + gene.chro)

                ##Correct coverage
                if gene.coverage != 1:
                    if gene.searchCDS(seq) is False:
                        sys.stderr.write("Gene " + gene.geneId() + " partial\n")
                        bad += 1
                        continue
                    else:
                        sys.stderr.write("Gene " + gene.geneId() + " correct\n")
                
                ##add sequence and MLST
                if gene.strand =="+":
                    sequence = seq.seq[gene.start:gene.end]
                else:
                    sequence = seq.seq[gene.start:gene.end].reverse_complement()
                    
                ##Verify sequence correct
                if len(sequence) != (gene.end-gene.start):
                    sys.stderr.write("Gene " + gene.geneId() + " incomplet\n")
                    bad += 1
                    continue
                try:
                    tmp = sequence.translate(table="Bacterial", cds=True)
                except:
                    sys.stderr.write("Gene " + gene.geneId() + " not correct ")
                    sys.stderr.write("Contig " + seq.id + " Start " + str(gene.start) + \
                                     " end " + str(gene.end) +"\n")
                    bad += 1
                    continue

                ##Insert data in database
                seqid = insert_sequence(cursor, str(sequence).upper())
                cursor2.execute('''INSERT INTO mlst(souche, gene, seqid) VALUES(?,?,?)''', \
                                    (name, gene.geneId(), seqid))
        db.commit()
        sys.stderr.write("Add " + str(len(genes)-bad) + " new MLST gene to database\n")
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
