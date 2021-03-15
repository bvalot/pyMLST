import logging
import os

from abc import ABC
from contextlib import contextmanager

import networkx as nx
from Bio import SeqIO

from pymlst.lib import blat, psl
from pymlst.wg_commands.db.database import DatabaseWG


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


def create_logger():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s: %(asctime)s] %(message)s')
    return logging.getLogger('wgMlst')


class WholeGenomeMLST:

    def __init__(self, file=None, ref='ref'):
        self.database = DatabaseWG(file)
        self.ref = ref
        self.blat_path = '/usr/bin/'  # TODO: change this
        self.logger = create_logger()

    def create(self, coregene, concatenate=True, remove=True):
        genes = set()
        to_remove = set()

        for gene in SeqIO.parse(coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID: " + gene.id)
            else:
                genes.add(gene.id)

            added, seq_id = self.database.add_sequence(str(gene.seq))

            if not added:
                if concatenate:
                    self.database.concatenate_gene(seq_id, gene.id)
                    self.logger.info("Concatenate gene " + gene.id + "\n")
                elif remove:
                    to_remove.add(seq_id)
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")
            else:
                self.database.add_mlst(self.ref, gene.id, seq_id)

        if to_remove:
            self.database.remove_sequences(to_remove)
            self.logger.info("Remove duplicate sequence: " + str(len(to_remove)) + "\n")

    def add_strain(self, strain, identity, coverage, genome):
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
            self.logger.info("Search coregene with BLAT\n")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage)
            self.logger.info("Finish run BLAT, found " + str(len(genes)) + " genes\n")

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
                            self.logger.info("Gene " + gene.geneId() + " partial: removed\n")
                            bad += 1
                            continue
                        else:
                            self.logger.info("Gene " + gene.geneId() + " fill: added\n")

                    # Verify CDS
                    if psl.testCDS(gene.getSequence(seq), False) is False:
                        if gene.searchCorrectCDS(seq, coverage) is False:
                            self.logger.info("Gene " + gene.geneId() + " not correct: removed\n")
                            bad += 1
                            continue
                        else:
                            self.logger.info("Gene " + gene.geneId() + " correct: added\n")

                    # add sequence and MLST
                    sequence = gene.getSequence(seq)

                    # Insert data in database
                    seqid = self.database.add_sequence(str(sequence))[1]
                    self.database.add_mlst(name, gene.geneId(), seqid)

            self.logger.info("Add " + str(len(genes) - bad) + " new MLST gene to database\n")
            self.logger.info("FINISH\n")

        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def remove_gene(self, list, genes):
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
                self.logger.info("Not found\n")
            else:
                self.logger.info("OK\n")

            self.database.remove_gene(gene)
            self.database.remove_orphan_sequences(seqids)

    def remove_strain(self, list, strains):
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
                self.logger.info("Not found\n")
            else:
                self.logger.info("OK\n")

            self.database.remove_strain(strain)
            self.database.remove_orphan_sequences(seqids)

    def extract(self, extractor, output):
        extractor.extract(self.database, self.ref, output, self.logger)

    def find(self, finder):
        finder.find(self.database)

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
    logger.info("Number of genes to look at : " + str(len(genes)) + "\n")

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
