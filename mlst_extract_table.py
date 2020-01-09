#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Extract MLST table from an wgMLST database"""

import sys
import os
import argparse
import sqlite3
from lib import __version__
import lib.sql as sql

desc = "Extract MLST table from a wgMLST database"
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
    help='Minimun number of strain found to keep a gene (default:0)')
command.add_argument('-k', '--keep', action='store_true', \
    help='Keep only gene with different allele (omit missing)')
command.add_argument('-d', '--duplicate', action='store_false', \
    help='Conserve duplicate gene (default remove)')
command.add_argument('-V', '--inverse', action='store_true', \
    help='Keep only gene that do not meet the filter of mincover or keep options')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to store MLST')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)
 
def get_all_strain(cursor):
    cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (sql.ref,))
    return [a[0] for a in cursor.fetchall()]

def get_all_gene(cursor):
    cursor.execute('''SELECT distinct(gene) FROM mlst WHERE souche = ?''', (sql.ref,))
    return [row[0] for row in cursor.fetchall()]

def get_duplicate_gene(cursor):
    cursor.execute('''SELECT gene
                      FROM mlst m
                      WHERE exists (
                      select 1 from mlst
                      where souche = m.souche
                      and gene = m.gene
                      and id != m.id
                      )
                      AND m.souche != ?
                      GROUP BY gene''', (sql.ref,))
    return [row[0] for row in cursor.fetchall()]

def get_count_seqid_by_gene(cursor):
    cursor.execute('''SELECT gene, count(distinct seqid)
                      from mlst
                      WHERE souche != ?
                      GROUP BY gene''', (sql.ref,))
    return {row[0]:row[1] for row in cursor.fetchall()}

def get_count_souche_by_gene(cursor):
    cursor.execute('''SELECT gene, count(distinct souche)
                      from mlst
                      WHERE souche != ?
                      GROUP by gene''', (sql.ref,))
    return {row[0]:row[1] for row in cursor.fetchall()}

def get_number_sequences(cursor):
    cursor.execute('''SELECT DISTINCT seqid FROM mlst WHERE souche!=?''', (sql.ref,))
    return len(cursor.fetchall())

def get_distance_between_strain(cursor, valid_shema):
    cursor.execute('''SELECT m1.souche, m2.souche, count( distinct m1.gene)
                      FROM mlst as m1 
                      JOIN mlst as m2 
                      ON m1.gene = m2.gene 
                      AND m1.souche != m2.souche 
                      AND m1.seqid != m2.seqid
                      WHERE m1.gene IN ( {} )
                      AND m1.souche != ?
                      AND m2.souche != ?
                      GROUP BY m1.souche, m2.souche'''.format(", ".join(["'" + g + "'" for g in valid_shema])), \
                   (sql.ref, sql.ref))
    distance = {}
    for row in cursor.fetchall():
        x = distance.setdefault(row[0], {})
        x[row[1]] = row[2]
    return distance

def get_mlst(cursor, valid_shema):
    cursor.execute('''SELECT gene, souche, group_concat(seqid) as seqid
                      FROM mlst
                      WHERE souche != ?
                      AND gene IN ( {} )
                      GROUP BY gene, souche'''.format(", ".join(["'" + g + "'" for g in valid_shema])), (sql.ref,))
    mlst = {}
    for row in cursor.fetchall():
        x = mlst.setdefault(row[0], {})
        x[row[1]] = row[2]
    return mlst

if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    output = args.output
    database = args.database
    export = args.export
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()

        ##index
        sql.index_database(cursor)

        ##read samples mlst
        strains = get_all_strain(cursor)
        ##Minimun number of strain
        if args.mincover < 0 or args.mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : " + str(len(strains)))

        ## allgene
        allgene = get_all_gene(cursor)
        
        ## duplicate gene
        dupli = set(get_duplicate_gene(cursor))

        ## cover without duplication
        count = get_count_souche_by_gene(cursor)

        ## Count distinct gene
        diff = get_count_seqid_by_gene(cursor)

        ##filter coregene that is not sufficient mincover or keep only different or return inverse
        valid_shema = []
        
        ## Test different case for validation
        for g in allgene:  
            valid = []
            if args.keep is True:
                if diff.get(g, 0) > 1:
                    valid.append(True)
                else:
                    valid.append(False)
            else:
                valid.append(True)
            if count.get(g, 0) >= args.mincover:
                valid.append(True)
            else:
                valid.append(False)
            if args.duplicate:
                if g in dupli:
                    valid.append(False)
                else:
                 valid.append(True)   
            else:
                valid.append(True)
            if args.inverse is False:
                if sum(valid) == 3:
                    valid_shema.append(g)
            else:
                if sum(valid) < 3:
                    valid_shema.append(g)

        ##report
        sys.stderr.write("Number of coregene used : " + str(len(valid_shema)) + \
                         "/" + str(len(allgene)) + "\n")
        
        ##export different case with choices
        if export == "strain":
            if args.count is False:
                output.write("\n".join(strains) +"\n")
            else:
                cursor.execute('''SELECT souche, count(distinct gene)
                                  FROM mlst
                                  WHERE gene IN ( {} )
                                  GROUP BY souche'''.format(", ".join(["'" + g + "'" for g in valid_shema])))
                tmp = {row[0]:row[1] for row in cursor.fetchall()}
                for strain in strains:
                    output.write(strain + "\t" + str(tmp.get(strain)) + "\n")
        elif export == "gene":
            output.write("\n".join(sorted(valid_shema)) + "\n")
        elif export == "distance":
            if args.duplicate is False:
                sys.stderr.write("WARNINGS : Calculate distance between strains " + \
                           "using duplicate genes could reported bad result\n")
            output.write(str(len(strains)) + "\n")
            distance = get_distance_between_strain(cursor, valid_shema)
            for s1 in strains:
                output.write(s1 + "\t")
                c = [str(distance.get(s1, {}).get(s2, 0)) for s2 in strains]
                output.write("\t".join(c) + "\n")
        elif export == "mlst":
            output.write("GeneId\t" + "\t".join(strains)+"\n")
            mlst = get_mlst(cursor, valid_shema)
            for g in valid_shema:
                towrite = [g]
                mlstg= mlst.get(g, {})
                for s in strains:
                    towrite.append(mlstg.get(s, ""))
                output.write("\t".join(towrite) + "\n")
        elif export == "stat":
            output.write("Strains\t"+str(len(strains))+"\n")
            output.write("Coregenes\t"+str(len(allgene))+"\n")
            output.write("Sequences\t"+str(get_number_sequences(cursor))+"\n")
        else:
            raise Exception("This export format is not supported: " + export)
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
