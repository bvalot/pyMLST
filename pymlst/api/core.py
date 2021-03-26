import logging
import os
import sys
import re

from abc import ABC
from contextlib import contextmanager

import networkx as nx
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError

from pymlst.lib.database import DatabaseCLA, DatabaseWG
from pymlst.lib import blat, psl


@contextmanager
def open_wg(file=None, ref='ref'):
    mlst = WholeGenomeMLST(file, ref)
    try:
        yield mlst
    except Exception:
        mlst.rollback()
        raise
    else:
        mlst.commit()
    finally:
        mlst.close()


@contextmanager
def open_cla(file=None, ref='ref'):
    mlst = ClassicalMLST(file, ref)
    try:
        yield mlst
    except Exception:
        mlst.rollback()
        raise
    else:
        mlst.commit()
    finally:
        mlst.close()


def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs


def strip_file(file):
    found = []
    if file is not None:
        for line in file.readLines():
            found.append(line.rstrip('\n'))
    return []


def compar_seqs(seqs):
    count = 0
    dim = len(seqs[0])
    for j in range(0, len(seqs[0])):
        d = set([s[j] for s in seqs])
        if '-' in d:
            d.remove('-')
        if len(d) > 1:
            count += 1
    return count


def write_count(count, texte):
    if count:
        count.write(texte)


def validate_sequence(sequence):
    try:
        sequence.translate(cds=True, table=11)
    except TranslationError:
        return False
    else:
        return True


def create_logger():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s: %(asctime)s] %(message)s')
    return logging.getLogger('mlst')


class ClassicalMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseCLA(file)
        self.ref = ref
        self.blat_path = '/usr/bin/'
        self.logger = create_logger()

    def create(self, scheme, alleles):
        # Verify sheme list with fasta files
        header = scheme.readline().rstrip("\n").split("\t")
        if len(header) != len(alleles) + 1:
            raise Exception("The number of genes in sheme don't correspond to the number of fasta file\n"
                            + " ".join(header) + "\n")
        fastas = {}
        for f in alleles:
            name = f.name.split("/")[-1]
            name = name[:name.rfind(".")]
            if name not in header:
                raise Exception("Gene " + name + " not found in sheme\n" + " ".join(header))
            fastas[name] = f

        # load sequence allele
        alleles = {}
        for g, f in fastas.items():
            alleles[g] = set()
            for seq in SeqIO.parse(f, 'fasta'):
                try:
                    # if len(seq.id.split("_")) == 2:
                    #     allele = int(seq.id.split("_")[1])
                    # elif len(seq.id.split("-")) == 2:
                    #     allele = int(seq.id.split("-")[1])
                    # elif g in seq.id:
                    #     allele = int(seq.id.replace(g, ""))
                    # else:
                    #     allele = int(seq.id)
                    match = re.search('[0-9]+$', seq.id)
                    allele = int(match.group(0))

                except Exception:
                    raise Exception("Unable to obtain allele number for the sequence: " + seq.id)
                self.database.add_sequence(str(seq.seq).upper(), g, allele)
                alleles.get(g).add(allele)

        # load MLST sheme
        for line in scheme:
            h = line.rstrip("\n").split("\t")
            st = int(h[0])
            for g, a in zip(header[1:], h[1:]):
                if int(a) not in alleles.get(g):
                    self.logger.info(
                        "Unable to find the allele number " + a + " for gene " + g + "; replace by 0")
                    self.database.add_mlst(st, g, 0)
                else:
                    self.database.add_mlst(st, g, int(a))

        self.logger.info('Database initialized')

    def search_st(self, genome, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 to 1")
        path = blat.test_blat_exe(self.blat_path)
        tmpfile, tmpout = blat.blat_tmp()

        try:
            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            self.logger.info("Search coregene with BLAT")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage, self.logger)
            self.logger.info("Finish run BLAT, found " + str(len(genes)) + " genes")

            # Search sequence MLST
            seqs = read_genome(genome)
            self.logger.info("Search allele gene to database")
            # print(genes)
            allele = {i: [] for i in coregenes}
            st = {i: set() for i in coregenes}
            for coregene in coregenes:
                if coregene not in genes:
                    allele.get(coregene).append("")
                    continue
                for gene in genes.get(coregene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # verify coverage and correct
                    if gene.coverage != 1:
                        gene.searchCorrect()
                        self.logger.info("Gene " + gene.geneId() + " fill: added")

                    # get sequence
                    sequence = str(gene.getSequence(seq)).upper()

                    # verify complet sequence
                    if len(sequence) != (gene.end - gene.start):
                        self.logger.info("Gene " + gene.geneId() + " removed")
                        continue

                    # write fasta file with coregene
                    if fasta is not None:
                        fasta.write(">" + coregene + "\n")
                        fasta.write(sequence + "\n")

                    # search allele
                    res = self.database.get_allele_by_sequence_and_gene(sequence, coregene)
                    if res is not None:
                        allele.get(coregene).append(str(res[0]))
                        # cursor.execute('''SELECT st FROM mlst WHERE gene=? and allele=?''',
                        #                (coregene, row[0]))
                        strains = self.database.get_strains_by_gene_and_allele(coregene, res[0])
                        for strain in strains:
                            st.get(coregene).add(strain[0])
                    else:
                        allele.get(gene.geneId()).append("new")

            # if only know allele or not found
            # Seach st
            st_val = []
            if sum([len(i) == 1 and i[0] != "new" for i in allele.values()]) == len(allele):
                tmp = None
                for s in st.values():
                    if s:
                        if tmp is None:
                            tmp = s
                        else:
                            tmp = tmp.intersection(s)
                st_val = list(tmp)

            # print result
            coregenes.sort()
            output.write("Sample\tST\t" + "\t".join(coregenes) + "\n")
            output.write(genome.name + "\t" + ";".join(map(str, st_val)))
            for coregene in coregenes:
                output.write("\t" + ";".join(map(str, allele.get(coregene))))
            output.write("\n")
            self.logger.info("FINISH")
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def __create_coregene(self, tmpfile):
        ref = int(1)
        all_rows = self.database.get_genes_by_allele(ref)
        coregenes = []
        for row in all_rows:
            tmpfile.write('>' + row[0] + "\n" + row[1] + "\n")
            coregenes.append(row[0])
        return coregenes

    def close(self):
        self.database.close()

    def commit(self):
        self.database.commit()

    def rollback(self):
        self.database.rollback()


class WholeGenomeMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseWG(file)
        self.ref = ref
        self.blat_path = '/usr/bin/'  # TODO: change this
        self.logger = create_logger()

    def create(self, coregene, concatenate=False, remove=False):
        genes = set()
        to_remove = set()
        rc_genes = 0
        invalid_genes = 0

        uh = 0
        for gene in SeqIO.parse(coregene, 'fasta'): # only 2503 elements in it
            uh += 1
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID: " + gene.id)
            else:
                genes.add(gene.id)

            if not validate_sequence(gene.seq):
                gene.seq = gene.seq.reverse_complement()
                if validate_sequence(gene.seq):
                    rc_genes += 1
                else:
                    invalid_genes += 1
                    continue

            added, seq_id = self.database.add_sequence(str(gene.seq))

            if not added:
                if concatenate:
                    self.database.concatenate_gene(seq_id, gene.id)
                    self.logger.info("Concatenate gene " + gene.id)
                elif remove:
                    to_remove.add(seq_id)
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")
            else:
                self.database.add_mlst(self.ref, gene.id, seq_id)

        print('Count: ' + str(uh))

        if to_remove:
            self.database.remove_sequences(to_remove)
            self.logger.info("Remove duplicate sequence: " + str(len(to_remove)))

        if rc_genes:
            self.logger.info('Reverse-complemented genes: ' + str(rc_genes))

        if invalid_genes:
            self.logger.info('Skipped invalid genes: ' + str(invalid_genes))

        self.logger.info('Database initialized')

    def add_strain(self, genome, strain=None, identity=0.95, coverage=0.90):
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 and 1")
        path = blat.test_blat_exe(self.blat_path)

        name = strain
        if name is None:
            name = genome.name.split('/')[-1]
        if ";" in name:
            raise Exception("Strain name must not contains special ';'\n")

        tmpfile, tmpout = blat.blat_tmp()

        try:
            # verify that the strain is not already in the database
            if self.database.contains_souche(name):
                raise Exception("Strain is already present in database:\n" + name)

            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            self.logger.info("Search coregene with BLAT")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage, self.logger)
            self.logger.info("Finish run BLAT, found " + str(len(genes)) + " genes")

            # add MLST sequence
            seqs = read_genome(genome)
            bad = 0
            for coregene in coregenes:
                if coregene not in genes:
                    continue
                for gene in genes.get(coregene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # Correct coverage
                    if gene.coverage != 1:
                        if gene.searchPartialCDS(seq, coverage) is False:
                            self.logger.info("Gene " + gene.geneId() + " partial: removed")
                            bad += 1
                            continue
                        else:
                            self.logger.info("Gene " + gene.geneId() + " fill: added")

                    # Verify CDS
                    if psl.testCDS(gene.getSequence(seq), False) is False:
                        if gene.searchCorrectCDS(seq, coverage) is False:
                            self.logger.info("Gene " + gene.geneId() + " not correct: removed")
                            bad += 1
                            continue
                        else:
                            self.logger.info("Gene " + gene.geneId() + " correct: added")

                    # add sequence and MLST
                    sequence = gene.getSequence(seq)

                    # Insert data in database
                    seqid = self.database.add_sequence(str(sequence))[1]
                    self.database.add_mlst(name, gene.geneId(), seqid)

            self.logger.info("Add " + str(len(genes) - bad) + " new MLST gene to databas")
            self.logger.info("FINISH")

        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def remove_gene(self, genes, list=None):
        # list genes to remove
        all_genes = strip_file(list)
        if genes is not None:
            all_genes.extend(genes)
        if len(all_genes) == 0:
            raise Exception("No gene to remove found.\n")
        all_genes = set(all_genes)

        for gene in all_genes:
            self.logger.info(gene + " : ")

            seqids = self.database.get_gene_sequences_ids(gene)
            if len(seqids) == 0:
                self.logger.info("Not found")
            else:
                self.logger.info("OK")

            self.database.remove_gene(gene)
            self.database.remove_orphan_sequences(seqids)

    def remove_strain(self, strains, list=None):
        if self.ref in strains:
            raise Exception("Ref schema could not be remove from this database")

        # list strains to remove
        all_strains = strip_file(list)
        if strains is not None:
            all_strains.extend(strains)
        if len(all_strains) == 0:
            raise Exception("No strain to remove found.\n")
        all_strains = set(all_strains)

        for strain in all_strains:
            self.logger.info(strain + " : ")

            seqids = self.database.get_strain_sequences_ids(strain)
            if len(seqids) == 0:
                self.logger.info("Not found")
            else:
                self.logger.info("OK")

            self.database.remove_strain(strain)
            self.database.remove_orphan_sequences(seqids)

    def extract(self, extractor, output=sys.stdout):
        extractor.extract(self.database, self.ref, output, self.logger)

    def __create_coregene(self, tmpfile):
        ref_genes = self.database.get_gene_by_souche(self.ref)
        coregenes = []
        for row in ref_genes:
            tmpfile.write('>' + row.gene + "\n" + row.sequence + "\n")
            coregenes.append(row[0])
        return coregenes

    def close(self):
        self.database.close()

    def commit(self):
        self.database.commit()

    def rollback(self):
        self.database.rollback()


class Extractor(ABC):
    def extract(self, base, ref, output, logger):
        pass


def find_recombination(genes, alignment, output):
    logger = create_logger()

    genes = [line.rstrip("\n") for line in genes]
    logger.info("Number of genes to look at : " + str(len(genes)))

    sequences = [[] for _ in genes]
    samples = []

    # load sequences by gene
    indice = 0
    for line in alignment:
        line = line.rstrip("\n")

        # header
        if line.startswith(">"):
            indice = 0
            samples.append(line.lstrip(">"))
            continue

        # check genes number correct
        if indice >= len(genes):
            raise Exception("The genes list seems not correspond to the alignment\n" + str(indice))

        # genes
        sequences[indice].append(line)
        indice += 1

    # check sequences are correctly align
    for i, seqs in enumerate(sequences):
        if len(set([len(s) for s in seqs])) > 1:
            print(set([len(s) for s in seqs]))
            raise Exception("Following genes seems to be not align: " + genes[i])

    for i, seqs in enumerate(sequences):
        c = compar_seqs(seqs)
        output.write(genes[i] + "\t" + str(c) + "\t" + str(len(seqs[0])) + "\n")


def find_subgraph(threshold, count, distance, output):
    samps = []
    dists = []
    try:
        strains = int(distance.readline().rstrip("\n"))
    except:
        raise Exception("The distance file seems not correctly formatted\n Not integer on first line")

    for line in distance.readlines():
        h = line.rstrip("\n").split("\t")
        samps.append(h[0])
        dists.append(h[1:])

    if len(samps) != strains:
        raise Exception("The distance file seems not correctly formatted\n Number of strains " + str(
            len(samps)) + " doesn't correspond to " + str(strains))

    # create graph
    G = nx.Graph()
    G.add_nodes_from(samps)

    for i, s in enumerate(samps):
        for j, d in enumerate(dists[i]):
            d = int(d)
            if i == j or d > threshold:
                continue
            G.add_edge(samps[i], samps[j], weight=d)

    # extract interconnected subgraph
    # count sample not found
    samps2 = set(samps)
    grps = []
    for subG in [G.subgraph(c) for c in nx.connected_components(G)]:

        inds = []
        for n in subG.nodes():
            samps2.remove(n)
            inds.append(samps.index(n))
        grps.append(inds)

    grps.sort(key=len, reverse=True)

    # write result
    write_count(count, "Group\t" + "\t".join(samps) + "\n")
    for i, g in enumerate(grps):
        a = len(samps) * [0]
        output.write("Group" + str(i))
        for n in g:
            a[n] = 1
            output.write(" " + samps[n])
        write_count(count, str(i) + "\t" + "\t".join(map(str, a)) + "\n")
        output.write("\n")

    if count:
        count.close()
    output.close()
