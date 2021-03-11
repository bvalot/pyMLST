#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Extract MLST table from an wgMLST database"""
import os
import sys

import click

from pymlst.lib import sql
import pandas as pd
import numpy as np
from pymlst.lib.benchmark import benchmark

from pymlst.wg_commands.db.database import DatabaseWG

desc = "Extract MLST table from a wgMLST database"

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
@benchmark
def cli(output, export, count, mincover, keep, duplicate, inverse, database):
    """Extract MLST table from a wgMLST database"""

    try:
        # db = sqlite3.connect(database.name)
        # cursor = db.cursor()
        #
        # # index
        # sql.index_database(cursor)
        db = DatabaseWG(os.path.abspath(database.name))

        # read samples mlst
        strains = db.get_all_strains(sql.ref)
        # Minimun number of strain
        if mincover < 0 or mincover > len(strains):
            raise Exception("Mincover must be between 0 to number of strains : " + str(len(strains)))

        # allgene
        allgene = db.get_all_genes(sql.ref)

        # duplicate gene
        dupli = db.get_duplicated_genes(sql.ref)

        # cover without duplication
        count_souches = db.count_souches_per_gene(sql.ref)

        # Count distinct gene
        diff = db.count_sequences_per_gene(sql.ref)

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
            if count_souches.get(g, 0) >= mincover:
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
                # cursor.execute('''SELECT souche, count(distinct gene)
                #                   FROM mlst
                #                   WHERE gene IN ( {} )
                #                   GROUP BY souche'''.format(", ".join(["'" + g + "'" for g in valid_shema])))
                # tmp = {row[0]:row[1] for row in cursor.fetchall()}
                tmp = db.count_genes_per_souche(valid_shema)
                for strain in strains:
                    output.write(strain + "\t" + str(tmp.get(strain)) + "\n")
        elif export == "gene":
            output.write("\n".join(sorted(valid_shema)) + "\n")
        elif export == "distance":
            if duplicate is False:
                sys.stderr.write("WARNINGS : Calculate distance between strains " +
                                 "using duplicate genes could reported bad result\n")
            output.write(str(len(strains)) + "\n")
            distance = db.get_strains_distances(sql.ref, valid_shema)
            for s1 in strains:
                output.write(s1 + "\t")
                c = [str(distance.get(s1, {}).get(s2, 0)) for s2 in strains]
                output.write("\t".join(c) + "\n")
        elif export == "mlst":
            output.write("GeneId\t" + "\t".join(strains)+"\n")
            mlst = db.get_mlst(sql.ref, valid_shema)
            for g in valid_shema:
                towrite = [g]
                mlstg= mlst.get(g, {})
                for s in strains:
                    towrite.append(mlstg.get(s, ""))
                output.write("\t".join(towrite) + "\n")
        elif export == "grapetree":
            mlst = db.get_mlst(sql.ref, valid_shema)
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
            output.write("Sequences\t"+str(db.get_sequences_number(sql.ref))+"\n")
        else:
            raise Exception("This export format is not supported: " + export)
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
