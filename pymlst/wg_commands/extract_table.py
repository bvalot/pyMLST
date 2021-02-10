#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Extract MLST table from an wgMLST database"""

import sys

import sqlite3
import click

from pymlst.lib import sql
import pandas as pd
import numpy as np

desc = "Extract MLST table from a wgMLST database"


def get_all_strain(cursor):
    cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche!=?''', (sql.ref,))
    return [a[0] for a in cursor.fetchall()]


def get_all_gene(cursor):
    cursor.execute('''SELECT distinct(gene) FROM mlst WHERE souche = ?''', (sql.ref,))
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
    return {row[0]: row[1] for row in cursor.fetchall()}


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
    cursor.execute('''SELECT gene, souche, group_concat(seqid, ";") as seqid
                      FROM mlst
                      WHERE souche != ?
                      AND gene IN ( {} )
                      GROUP BY gene, souche'''.format(", ".join(["'" + g + "'" for g in valid_shema])), (sql.ref,))
    mlst = {}
    for row in cursor.fetchall():
        x = mlst.setdefault(row[0], {})
        x[row[1]] = row[2]
    return mlst


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Export MLST table to (default=stdout)')
@click.option('--export', '-e',
              type=click.Choice(['mlst', 'grapetree', 'distance', 'strain', 'gene', 'stat']),
              default='mlst',
              help='Defined the export format')
@click.option('--count', '-c',
              is_flag=True,
              help='In strain mode, count the number of gene present in the database')
@click.option('--mincover', '-m',
              type=int, default=0,
              help='Minimun number of strain found to keep a gene (default:0)')
@click.option('--keep', '-k',
              is_flag=True,
              help='Keep only gene with different allele (omit missing)')
@click.option('--duplicate', '-d',
              is_flag=True, default=True,
              help='Conserve duplicate gene (default remove)')
@click.option('--inverse', '-V',
              is_flag=True,
              help='Keep only gene that do not ' \
              'meet the filter of mincover or keep options')
@click.argument('database', type=click.File('r'))
def cli(output, export, count, mincover, keep, duplicate, inverse, database):
    """Extract MLST table from a wgMLST database"""
    print('dupli: ', duplicate)

    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()

        # index
        sql.index_database(cursor)

        # read samples mlst
        strains = get_all_strain(cursor)
        # Minimun number of strain
        if mincover < 0 or mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : " + str(len(strains)))

        # allgene
        allgene = get_all_gene(cursor)

        # duplicate gene
        dupli = sql.get_duplicate_gene(cursor)

        # cover without duplication
        count = get_count_souche_by_gene(cursor)

        # Count distinct gene
        diff = get_count_seqid_by_gene(cursor)

        # filter coregene that is not sufficient mincover or keep only different or return inverse
        valid_shema = []

        # Test different case for validation
        for g in allgene:
            valid = []
            if keep is True:
                if diff.get(g, 0) > 1:
                    valid.append(True)
                else:
                    valid.append(False)
            else:
                valid.append(True)
            if count.get(g, 0) >= mincover:
                valid.append(True)
            else:
                valid.append(False)
            if duplicate:
                if g in dupli:
                    valid.append(False)
                else:
                    valid.append(True)
            else:
                valid.append(True)
            if inverse is False:
                if sum(valid) == 3:
                    valid_shema.append(g)
            else:
                if sum(valid) < 3:
                    valid_shema.append(g)

        # report
        sys.stderr.write("Number of coregene used : " + str(len(valid_shema)) + \
                         "/" + str(len(allgene)) + "\n")

        # export different case with choices
        if export == "strain":
            if count is False:
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
            if duplicate is False:
                sys.stderr.write("WARNINGS : Calculate distance between strains " +
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
        elif export == "grapetree":
            mlst = get_mlst(cursor, valid_shema)
            df = pd.DataFrame(columns=["#GeneId"] + strains)
            for g in valid_shema:
                row = {"#GeneId": g}
                mlstg= mlst.get(g, {})
                for s in strains:
                    row[s] = mlstg.get(s, np.NaN)
                df = df.append(row, ignore_index=True)
            df = df.set_index('#GeneId')
            df = df.transpose()
            df = df.fillna(-1).astype(int)
            df.to_csv(output, sep='\t')
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
