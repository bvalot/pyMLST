"""Core classes and functions to work with Whole Genome MLST data."""

import logging
import os
import sys
from abc import ABC

import networkx as nx
from Bio import SeqIO
from decorator import contextmanager
from diff_match_patch import diff_match_patch
from sqlalchemy.sql.functions import count

from sqlalchemy import create_engine, bindparam, not_, Index
from sqlalchemy import and_
from sqlalchemy import func
from sqlalchemy import MetaData, Table, Column, Integer, Text, ForeignKey

from sqlalchemy.sql import select, exists
from sqlalchemy.sql import distinct
from sqlalchemy.sql.operators import in_op as in_

from pymlst.common import blat, psl
from pymlst.common.binaries import get_binary_path
from pymlst.common import utils


@contextmanager
def open_wg(file=None, ref='ref'):
    """A context manager function to wrap the creation a
       :class:`~pymlst.wg.core.WholeGenomeMLST` object.

    Context managers allow you to instantiate objects using the ``with`` keyword,
    this way you don't have to manage exceptions and the committing/closing processes yourself.

    :param file: The path to the database file to work with.
    :param ref: The name that will be given to the reference strain in the database.

    :yields: A :class:`~pymlst.wg.core.WholeGenomeMLST` object.
    """
    mlst = WholeGenomeMLST(file, ref)
    try:
        yield mlst
    except Exception:
        mlst.rollback()
        raise
    else:
        mlst.commit(renew=False)
    finally:
        mlst.close()


class DatabaseWG:
    """A core level class to manipulate the genomic database.

        .. warning:: Shouldn't be instantiated directly,
                     see :class:`~pymlst.wg.core.WholeGenomeMLST` instead.
    """

    def __init__(self, path=None):
        """
        :param path: The path to the database file to work with.
        """
        if path is None:
            self.engine = create_engine('sqlite://')  # creates a :memory: database
        else:
            self.engine = create_engine('sqlite:///' + path)  # must be an absolute path

        metadata = MetaData()

        self.sequences = Table('sequences', metadata,
                               Column('id', Integer, primary_key=True),
                               Column('sequence', Text),
                               Column('gene', Text))

        Index('seq_gene_ind', self.sequences.c.sequence, self.sequences.c.gene, unique=True)

        self.mlst = Table('mlst', metadata,
                          Column('id', Integer, primary_key=True),
                          Column('souche', Text),
                          Column('gene', Text),
                          Column('seqid', Integer, index=True))

        Index('souche_gene_ind', self.mlst.c.souche, self.mlst.c.gene, unique=True)

        metadata.create_all(self.engine)

        self.connection = self.engine.connect()
        self.transaction = self.connection.begin()

        self.cached_queries = {}

    def __get_cached_query(self, name, query_supplier):
        if name in self.cached_queries:
            return self.cached_queries[name]
        query = query_supplier()
        self.cached_queries[name] = query
        return query

    def get_sequences_by_gene(self, gene):
        return self.connection.execute(
            select([self.mlst.c.souche, self.sequences.c.sequence])
            .where(and_(
                self.mlst.c.seqid == self.sequences.c.id,
                self.mlst.c.gene == gene))
        ).fetchall()

    def add_mlst(self, souche, gene, seqid):
        """Adds an MLST gene bound to an existing sequence."""
        self.connection.execute(
            self.mlst.insert(),
            souche=souche, gene=gene, seqid=seqid)

    def add_sequence(self, sequence, gene):
        """Adds a sequence."""
        # query = self.__get_cached_query(
        #     'add_sequence',
        #     lambda:
        #     select([self.sequences.c.id])
        #         .where(and_(
        #             self.sequences.c.sequence == bindparam('sequence'),
        #             self.sequences.c.gene == bindparam('gene'))))
        #
        # existing = self.connection.execute(
        #     query,
        #     sequence=sequence,
        #     gene=gene
        # ).fetchone()
        #
        # if existing is not None:
        #     return False, existing.id
        res = self.connection.execute(
            self.sequences.insert(),
            sequence=sequence,
            gene=gene)
        # return True, res.inserted_primary_key[0]
        return res.inserted_primary_key[0]

    def concatenate_gene(self, seq_id, gene_name):
        """Associates a new gene to an existing sequence using concatenation."""
        self.connection.execute(
            self.mlst.update()
                .values(gene=self.mlst.c.gene + ';' + gene_name)
                .where(self.mlst.c.seqid == seq_id))

    def remove_sequences(self, ids):
        """Removes sequences and the genes that reference them."""
        self.connection.execute(
            self.sequences.delete()
            .where(in_(self.sequences.c.id, ids)))
        self.connection.execute(
            self.mlst.delete()
            .where(in_(self.mlst.c.seqid, ids)))

    def remove_orphan_sequences(self, ids):
        """Removes sequences if they aren't referenced by any gene."""
        query = self.sequences.delete() \
            .where(and_(
                not_(exists(
                    select([self.mlst.c.id])
                    .where(self.mlst.c.seqid == self.sequences.c.id))),
                self.sequences.c.id == bindparam('seqid')))

        for seqid in ids:
            self.connection.execute(
                query,
                seqid=seqid[0])

    def remove_gene(self, gene):
        """Removes a specific gene."""
        self.connection.execute(
            self.mlst.delete()
            .where(self.mlst.c.gene == gene))

    def remove_strain(self, strain):
        """Removes a specific strain."""
        self.connection.execute(
            self.mlst.delete()
                .where(self.mlst.c.souche == strain))

    def get_seq_id_by_sequence(self, sequence):
        """Gets a sequence ID, searching by sequence."""
        return self.connection.execute(
            select([self.sequences.c.id])
            .where(self.sequences.c.sequence == sequence)
        ).fetchone()

    def get_seq_id_by_sequence_and_gene(self, sequence, gene):
        """Gets a sequence ID, searching by sequence and gene."""
        return self.connection.execute(
            select([self.sequences.c.id])
            .where(and_(
                self.sequences.c.sequence == sequence,
                self.sequences.c.gene == gene))
        ).fetchone()

    def get_gene_sequences_ids(self, gene):
        """Gets the IDs of the sequences associated with a specific gene."""
        return self.connection.execute(
            select([self.mlst.c.seqid])
            .where(self.mlst.c.gene == gene)
        ).fetchall()

    def get_strain_sequences_ids(self, strain):
        """Gets the IDs of the sequences associated with a specific strain."""
        return self.connection.execute(
            select([self.mlst.c.seqid])
            .where(self.mlst.c.souche == strain)
        ).fetchall()

    def get_sequence_by_gene_and_strain(self, gene, souche):
        """Gets the sequence associated to a specific gene in a specific strain."""
        return self.connection.execute(
            select([self.mlst.c.gene, self.sequences.c.sequence, self.sequences.c.id])
            .where(and_(
                self.mlst.c.gene == gene,
                self.mlst.c.souche == souche,
                self.mlst.c.seqid == self.sequences.c.id))
        ).fetchone()

    def get_sequences_by_souche(self, souche):
        """Gets all the sequences associated to a specific strain."""
        return self.connection.execute(
            select([self.mlst.c.gene, self.sequences.c.sequence])
            .where(and_(
                self.mlst.c.souche == souche,
                self.mlst.c.seqid == self.sequences.c.id))
        ).fetchall()

    def contains_souche(self, souche):
        """Whether the strain exists in the base or not."""
        return self.connection.execute(
            select([self.mlst.c.id])
            .where(self.mlst.c.souche == souche)
            .limit(1)
        ).fetchone() is not None

    def get_gene_sequences(self, gene, ref):
        """Gets all the sequences for a specific gene and
        lists the strains that are referencing them."""
        query = self.__get_cached_query(
            'get_gene_sequences',
            lambda:
            select([self.mlst.c.seqid,
                    func.group_concat(self.mlst.c.souche, bindparam('separator')),
                    self.sequences.c.sequence])
            .select_from(
                self.mlst.join(self.sequences))
            .where(and_(
                self.mlst.c.souche != bindparam('ref'),
                self.mlst.c.gene == bindparam('gene')))
            .group_by(self.mlst.c.seqid))

        res = self.connection.execute(
            query,
            separator=';',
            ref=ref,
            gene=gene
        ).fetchall()

        seqs = []
        for seq in res:
            tmp = seq[1].split(";")
            tmp.sort()
            seqs.append([seq[0], tmp, seq[2]])
        return seqs

    def get_genes_coverages(self, ref):
        """Counts the number of strains referencing each gene."""
        return self.connection.execute(
            select([self.mlst.c.gene,
                    func.count(distinct(self.mlst.c.souche))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()

    def get_duplicated_genes(self, ref):
        """Gets the genes that are duplicated."""
        m_alias = self.mlst.alias()

        exist_sub = select([self.mlst]) \
            .where(and_(
                self.mlst.c.souche == m_alias.c.souche,
                self.mlst.c.gene == m_alias.c.gene,
                self.mlst.c.id != m_alias.c.id))

        res = self.connection.execute(
            select([self.mlst.c.gene])
            .where(and_(
                exists(exist_sub),
                self.mlst.c.souche != ref))
            .group_by(self.mlst.c.gene)
        ).fetchall()

        return {row[0] for row in res}

    def get_all_strains(self, ref):
        """Gets all distinct strains."""
        res = self.connection.execute(
            select([distinct(self.mlst.c.souche)]).
            where(self.mlst.c.souche != ref)
        ).fetchall()
        return [r[0] for r in res]

    def get_all_genes(self, ref):
        """Gets all distinct genes."""
        res = self.connection.execute(
            select([distinct(self.mlst.c.gene)]).
            where(self.mlst.c.souche == ref)
        ).fetchall()
        return [r[0] for r in res]

    def count_sequences_per_gene(self, ref):
        """Gets the number of distinct sequences per gene."""
        res = self.connection.execute(
            select([self.mlst.c.gene, count(distinct(self.mlst.c.seqid))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_souches_per_gene(self, ref):
        """Gets the number of distinct stains per gene."""
        res = self.connection.execute(
            select([self.mlst.c.gene, count(distinct(self.mlst.c.souche))])
            .where(self.mlst.c.souche != ref)
            .group_by(self.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_genes_per_souche(self, valid_shema):
        """Gets the number of distinct genes per strain."""
        res = self.connection.execute(
            select([self.mlst.c.souche, count(distinct(self.mlst.c.gene))])
            .where(in_(self.mlst.c.gene, valid_shema))
            .group_by(self.mlst.c.souche)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def get_sequences_number(self, ref):
        """Gets the number of distinct."""
        return self.connection.execute(
               select([count(distinct(self.mlst.c.seqid))])
               .where(self.mlst.c.souche != ref)
        ).fetchone()[0]

    def get_strains_distances(self, ref, valid_schema):
        """Gets the strains distances.

        For all the possible pairs of strains, counts how many of their genes
        are different (different seqids so different sequences)
        The compared genes are restricted to the ones given in the valid_schema.
        """
        alias_1 = self.mlst.alias()
        alias_2 = self.mlst.alias()

        result = self.connection.execute(
            select(
                [alias_1.c.souche, alias_2.c.souche, count(distinct(alias_1.c.gene))])
            .select_from(
                alias_1.join(
                    alias_2,
                    and_(
                        alias_1.c.gene == alias_2.c.gene,
                        alias_1.c.souche != alias_2.c.souche,
                        alias_1.c.seqid != alias_2.c.seqid)))
            .where(
                and_(
                    in_(alias_1.c.gene, valid_schema),
                    alias_1.c.souche != ref,
                    alias_2.c.souche != ref))
            .group_by(alias_1.c.souche, alias_2.c.souche)
        ).fetchall()

        distance = {}
        for entry in result:
            dist = distance.setdefault(entry[0], {})
            dist[entry[1]] = entry[2]

        return distance

    def get_mlst(self, ref, valid_schema):
        """Gets the MLST sequences and their strains associated to the genes in the given schema."""
        result = self.connection.execute(
            select([self.mlst.c.gene, self.mlst.c.souche,
                    func.group_concat(self.mlst.c.seqid, ';')])
            .where(and_(self.mlst.c.souche != ref,
                        in_(self.mlst.c.gene, valid_schema)))
            .group_by(self.mlst.c.gene, self.mlst.c.souche)
        ).fetchall()

        mlst = {}

        for entry in result:
            sequences = mlst.setdefault(entry[0], {})
            sequences[entry[1]] = entry[2]
        return mlst

    def commit(self, renew):
        """Commits the modifications."""
        self.transaction.commit()
        if renew:
            self.transaction = self.connection.begin()

    def rollback(self):
        """Rollback the modifications."""
        self.transaction.rollback()

    def close(self):
        """Closes the database engine."""
        self.engine.dispose()


class WholeGenomeMLST:
    """Whole Genome MLST python representation.

    Example of usage::

        open_wg('database.db') as db:
            db.create(open('genome.fasta'))
            db.add_strain(open('strain_1.fasta'))
            db.add_strain(open('strain_2.fasta'))
    """

    def __init__(self, file, ref):
        """
        :param file: The path to the database file to work with.
        :param ref: The name that will be given to the reference strain in the database.
        """
        self.database = DatabaseWG(file)
        self.ref = ref
        self.patcher = diff_match_patch()

        utils.create_logger()

    def __add_mlst(self, sequence, gene, strain):
        ref_gene = self.database.get_sequence_by_gene_and_strain(gene, self.ref)
        patches = self.patcher.patch_make(ref_gene.sequence, sequence)
        if len(patches) == 0:
            new_id = ref_gene.id
        else:
            diff = self.patcher.patch_toText(patches)
            patch_id = self.database.get_seq_id_by_sequence_and_gene(diff, gene)  # test if patch already exists
            if patch_id is None:
                new_id = self.database.add_sequence(diff, gene)
            else:
                new_id = patch_id[0]
        self.database.add_mlst(strain, gene, new_id)

    def create(self, coregene, concatenate=False, remove=False):
        """Creates a whole genome MLST database from a core genome `fasta`_ file.

        :param coregene: The fasta file containing the reference core genome.
        :param concatenate: Whether we should concatenate genes with identical sequences.
        :param remove: Whether we should remove genes with identical sequences.

        For instance, if concatenate is set to ``True``, 2 genes **g1** and **g2**
        having the same sequence
        will be stored as a single gene named **g1;g2**.
        """
        genes = set()
        to_remove = set()
        rc_genes = 0
        invalid_genes = 0

        for gene in SeqIO.parse(coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID: " + gene.id)
            genes.add(gene.id)

            if not utils.validate_sequence(gene.seq):
                gene.seq = gene.seq.reverse_complement()
                if utils.validate_sequence(gene.seq):
                    rc_genes += 1
                else:
                    invalid_genes += 1
                    continue

            #added, seq_id = self.database.add_sequence(str(gene.seq))
            sequence = str(gene.seq)
            seq_id = self.database.get_seq_id_by_sequence(sequence)

            if seq_id is None: # This is a new sequence
                seq_id = self.database.add_sequence(sequence, gene.id)
            else: # A similar sequence was found
                if concatenate:
                    self.database.concatenate_gene(seq_id, gene.id)
                    logging.info("Concatenate gene %s", gene.id)
                elif remove:
                    to_remove.add(seq_id)
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")

            self.database.add_mlst(self.ref, gene.id, seq_id)

        if to_remove:
            self.database.remove_sequences(to_remove)
            logging.info("Remove duplicate sequence: %s", str(len(to_remove)))

        if rc_genes:
            logging.info('Reverse-complemented genes: %s', str(rc_genes))

        if invalid_genes:
            logging.info('Skipped invalid genes: %s', str(invalid_genes))

        logging.info('Database initialized')

    def add_strain(self, genome, strain=None, identity=0.95, coverage=0.90):
        """Adds a strain to the database.

        How it works:

        1. A `BLAT`_ research is performed on each given contig of the strain to find
           sub-sequences matching the core genes.
        2. The identified sub-sequences are extracted and added to our database where
           they are associated to a **sequence ID**.
        3. An MLST entry is created, referencing the sequence,
           the gene it belongs to, and the strain it was found in.

        :param genome: The strain genome we want to add as a `fasta`_ file.
        :param strain: The name that will be given to the new strain in the database.
        :param identity: Sets the minimum identity used by `BLAT`_ for sequences research (in percent).
        :param coverage: Sets the minimum accepted coverage for found sequences.
        """
        if identity < 0 or identity > 1:
            raise Exception("Identity must be between 0 and 1")

        path = get_binary_path('blat')
        if path is None:
            raise Exception('Unable to locate the Blat executable\n')

        name = strain
        if name is None:
            name = genome.name.split('/')[-1]
        if ";" in name:
            raise Exception("Strain name must not contains special ';'\n")

        tmpfile, tmpout = blat.blat_tmp()
        tmpout.close()

        try:
            # verify that the strain is not already in the database
            if self.database.contains_souche(name):
                raise Exception("Strain is already present in database:\n" + name)

            # read coregene
            coregenes = self.__create_coregene(tmpfile)
            tmpfile.close()

            # BLAT analysis
            logging.info("Search coregene with BLAT")
            genes = blat.run_blat(path, genome, tmpfile, tmpout, identity, coverage)
            logging.info("Finish run BLAT, found %s genes", str(len(genes)))

            # add MLST sequence
            seqs = utils.read_genome(genome)

            bad = 0
            partial = 0
            partial_filled = 0

            for coregene in coregenes:

                if coregene not in genes:
                    continue

                for gene in genes.get(coregene):

                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # Correct coverage
                    if gene.coverage == 1:
                        sequence = gene.get_sequence(seq)
                    else:
                        coregene_seq = self.database.get_sequence_by_gene_and_strain(
                            coregene, self.ref)[1]
                        sequence = gene.get_aligned_sequence(seq, coregene_seq)

                    if sequence and psl.test_cds(sequence):
                        if gene.coverage != 1:
                            logging.debug("Gene %s fill: added", gene.gene_id())
                            partial_filled += 1
                            partial += 1
                    else:
                        if gene.coverage == 1:
                            logging.debug("Gene %s not correct: removed", gene.gene_id())
                        else:
                            logging.debug("Gene %s partial: removed", gene.gene_id())
                            partial += 1
                        bad += 1
                        continue

                    # Insert data in database
                    self.__add_mlst(str(sequence), gene.gene_id(), name)
                    # seqid = self.database.add_sequence(str(sequence))[1]
                    # self.database.add_mlst(name, gene.gene_id(), seqid)

            logging.info("Added %s new MLST genes to the database", len(genes) - bad)
            logging.info('Found %s partial genes, filled %s', partial, partial_filled)
            logging.info('Removed %s genes', bad)
            logging.info("DONE")

        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)

    def remove_gene(self, genes, file=None):
        """Removes genes from the database.

        :param genes: Names of the genes to remove.
        :param file: A file containing a gene name per line.
        """
        # list genes to remove
        all_genes = utils.strip_file(file)
        if genes is not None:
            all_genes.extend(genes)
        if len(all_genes) == 0:
            raise Exception("No gene to remove found.\n")
        all_genes = set(all_genes)

        for gene in all_genes:
            logging.info("%s : ", gene)

            seqids = self.database.get_gene_sequences_ids(gene)
            if len(seqids) == 0:
                logging.info("Not found")
            else:
                logging.info("OK")

            self.database.remove_gene(gene)
            self.database.remove_orphan_sequences(seqids)

    def remove_strain(self, strains, file=None):
        """Removes entire strains from the database.

        :param strains: Names of the strains to remove.
        :param file: A file containing a strain name per line.
        """
        if self.ref in strains:
            raise Exception("Ref schema could not be remove from this database")

        # list strains to remove
        all_strains = utils.strip_file(file)
        if strains is not None:
            all_strains.extend(strains)
        if len(all_strains) == 0:
            raise Exception("No strain to remove found.\n")
        all_strains = set(all_strains)

        for strain in all_strains:
            logging.info("%s : ", strain)

            seqids = self.database.get_strain_sequences_ids(strain)
            if len(seqids) == 0:
                logging.info("Not found")
            else:
                logging.info("OK")

            self.database.remove_strain(strain)
            self.database.remove_orphan_sequences(seqids)

    def extract(self, extractor, output=sys.stdout):
        """Takes an extractor object and writes the extraction result on the given output.

        :param extractor: A :class:`~pymlst.wg.core.Extractor` object describing
                          the way data should be extracted.
        :param output: The output that will receive extracted data.
        """
        extractor.extract(self.database, self.ref, output)

    def __create_coregene(self, tmpfile):
        ref_genes = self.database.get_sequences_by_souche(self.ref)
        coregenes = []
        for row in ref_genes:
            tmpfile.write('>' + row.gene + "\n" + row.sequence + "\n")
            coregenes.append(row[0])
        return coregenes

    def commit(self, renew=True):
        """Commits the modifications."""
        self.database.commit(renew)

    def rollback(self):
        """Rollback the modifications."""
        self.database.rollback()

    def close(self):
        """Closes the database engine."""
        self.database.close()


class Extractor(ABC):
    """A simple interface to ease the process of creating new extractors."""
    def extract(self, base, ref, output):
        """
        :param base: The database to extract data from.
        :param ref: The name of the reference genome.
        :param output: The output where to write the extraction results.
        """


def find_recombination(genes, alignment, output):
    """Counts the number of versions of each gene.

    :param genes: List of genes (output of :class:`~pymlst.wg.extractors.TableExtractor`
                  using ``export='gene'``).
    :param alignment: `fasta`_ file alignment
                      (output of :class:`~pymlst.wg.extractors.SequenceExtractor` using ``align=True``).
    :param output: The output where to write the results.
    """
    genes = [line.rstrip("\n") for line in genes]
    logging.info("Number of genes to look at : %s", len(genes))

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
        if len({len(s) for s in seqs}) > 1:
            print({len(s) for s in seqs})
            raise Exception("Following genes seems to be not align: " + genes[i])

    for i, seqs in enumerate(sequences):
        compared = utils.compar_seqs(seqs)
        output.write(genes[i] + "\t" + str(compared) + "\t" + str(len(seqs[0])) + "\n")


def find_subgraph(distance, threshold=50, output=sys.stdout, export='group'):
    """Searches groups of strains separated by a distance threshold.

    :param threshold: Minimum distance to maintain for groups extraction.
    :param distance: Distance matrix file
                     (output of :class:`~pymlst.wg.extractors.TableExtractor`
                     with ``export='distance'``).
    :param output: The output where to write the results.
    :param export: Sets the export type.
    """
    samps = []
    dists = []
    try:
        strains = int(distance.readline().rstrip("\n"))
    except Exception as err:
        raise Exception("The distance file seems not correctly formatted\n Not integer on first line") from err

    for line in distance.readlines():
        dist_line = line.rstrip("\n").split("\t")
        samps.append(dist_line[0])
        dists.append(dist_line[1:])

    if len(samps) != strains:
        raise Exception("The distance file seems not correctly formatted\n Number of strains " +
                        str(len(samps)) + " doesn't correspond to " + str(strains))

    # create graph
    graph = nx.Graph()
    graph.add_nodes_from(samps)

    for strain_index, _ in enumerate(samps):
        for dist_index, dist in enumerate(dists[strain_index]):
            dist = int(dist)
            if strain_index == dist_index or dist > threshold:
                continue
            graph.add_edge(samps[strain_index], samps[dist_index], weight=dist)

    # extract interconnected subgraph
    # count sample not found
    samps2 = set(samps)
    grps = []
    for sub_graph in [graph.subgraph(c) for c in nx.connected_components(graph)]:

        inds = []
        for node in sub_graph.nodes():
            samps2.remove(node)
            inds.append(samps.index(node))
        grps.append(inds)

    grps.sort(key=len, reverse=True)

    # write result
    if export == 'group':
        for i, group in enumerate(grps):
            output.write('Group' + str(i))
            for node in group:
                output.write(" " + samps[node])
            output.write("\n")

    elif export == 'count':
        output.write('Group\t' + '\t'.join(samps) + '\n')
        for i, group in enumerate(grps):
            line = len(samps) * [0]
            for node in group:
                line[node] = 1
            output.write(str(i) + '\t' + '\t'.join(map(str, line)) + '\n')
    else:
        for i, group in enumerate(grps):
            for node in group:
                output.write('Group' + str(i) + '\t' + samps[node] + '\n')
