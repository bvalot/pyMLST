#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Extract MLST table from an wgMLST database"""

import sys
import os
import argparse
import sqlite3

ref="ref"

desc = "Extract MLST table from an wgMLST database"
command = argparse.ArgumentParser(prog='mlst_extract_table.py', \
    description=desc, usage='%(prog)s [options] database')
command.add_argument('-o', '--output', nargs='?', \
    type=argparse.FileType("w"), default=sys.stdout, \
    help='Export MLST table to (default=stdout)')
command.add_argument('-e', '--export', default="mlst", \
    choices=("mlst", "distance", "strain", "gene", "stat"), \
    help="; ".join(["Defined the export format", "mlst: MLST table (Default)",\
                   "distance: the distance matrix", "strain: the strain list", \
                    "gene: the gene list", "stat: statistics of the database"]))
command.add_argument('-c', '--count', action='store_true', \
    help='In strain mode, count the number of gene present in the database')
command.add_argument('-m', '--mincover', nargs='?', \
    type=int, default=0, \
    help='Minimun number of strain found to conserved a gene (default:0)')
command.add_argument('-k', '--keep', action='store_true', \
    help='Keep only gene with different allele (omit missing)')
command.add_argument('-V', '--inverse', action='store_true', \
    help='Conserved only gene that not pass filter of mincover or keep options')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to stock MLST')

def get_mlst(cursor, strain, shema):
    coregenes = [[] for i in range(0,len(shema))]
    cursor.execute('''SELECT gene,seqid FROM mlst WHERE souche=?''', (strain,))
    for row in cursor.fetchall():
        coregenes[shema.index(row[0])].append(row[1])
    return coregenes

def get_shema(cursor):
    ref = "ref"
    cursor.execute('''SELECT gene FROM mlst WHERE souche=?''', (ref,))
    return [a[0] for a in cursor.fetchall()]
 
def get_all_strain(cursor):
    cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (ref,))
    return [a[0] for a in cursor.fetchall()]

def get_number_sequences(cursor):
    cursor.execute('''SELECT DISTINCT seqid FROM mlst WHERE souche!=?''', (ref,))
    return len(cursor.fetchall())
    
def compare_strain(strain1, strain2, valid_shema):
    count = 0
    for a,(i,j) in enumerate(zip(strain1, strain2)):
        if a not in valid_shema:
            continue
        if i and j:
            if len(i) == len(j):
                if i!=j:
                    count +=1
            else:
                if not set(i).intersection(set(j)):
                    count += 1
    return count
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    output = args.output
    database = args.database
    export = args.export
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()

        ##read samples mlst
        strains = get_all_strain(cursor)
        ##Minimun number of strain
        if args.mincover < 0 or args.mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : " + str(len(strains)))

        ##read shema
        shema = get_shema(cursor)

        mlst = {}
        for strain in strains:
            mlst[strain] = get_mlst(cursor, strain, shema)

        ##filter coregene that is not sufficient mincover or keep only different or return inverse
        valid_shema = []
        for i,g in enumerate(shema):
            count = 0
            tmp = set()
            for s in strains:
                val = mlst.get(s)[i]
                if val:
                    count += 1
                    tmp = tmp.union(set(val))
            ## Test different case for validation
            valid = []
            if args.keep is True:
                if len(tmp) > 1:
                    valid.append(True)
                else:
                    valid.append(False)
            else:
                valid.append(True)
            if count >= args.mincover:
                valid.append(True)
            else:
                valid.append(False)
            if args.inverse is False:
                if sum(valid) == 2:
                    valid_shema.append(i)
            else:
                if sum(valid) < 2:
                    valid_shema.append(i)                
        
        ##export different case with choices
        if export == "strain":
            if args.count is False:
                output.write("\n".join(strains) +"\n")
            else:
                for strain in strains:
                    output.write(strain+"\t")
                    count = 0
                    for i in valid_shema:
                        val = mlst.get(strain)[i]
                        if val:
                            count += 1
                    output.write(str(count)+"\n")  
        elif export == "gene": 
            output.write("\n".join([g for i,g in enumerate(shema) if i in valid_shema]) + "\n")
        elif export == "distance":
            #output.write("\t" + "\t".join(strains)+"\n")
            output.write(str(len(strains)) + "\n")
            for s1 in strains:
                output.write(s1 + "\t")
                c = [str(compare_strain(mlst.get(s1), mlst.get(s2), valid_shema)) for s2 in strains]
                output.write("\t".join(c) + "\n")
        elif export == "mlst":
            output.write("GeneId\t" + "\t".join(strains)+"\n")
            for i,g in enumerate(shema):
                if i not in valid_shema:
                    continue
                towrite = [g]
                for s in strains:
                    val = mlst.get(s)[i]
                    if not val:
                        towrite.append("")
                    else:
                        towrite.append(";".join(map(str,val)))
                output.write("\t".join(towrite) + "\n")
        elif export == "stat":
            output.write("Strains\t"+str(len(strains))+"\n")
            output.write("Coregenes\t"+str(len(shema))+"\n")
            output.write("Sequences\t"+str(get_number_sequences(cursor))+"\n")
        else:
            raise Exception("This export format is not supported: " + export)
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
