""""Core classes and functions to work with Whole Genome MLST data."""
import logging
import os
import sys
from abc import ABC, abstractmethod
from enum import Enum

import networkx as nx
from Bio import SeqIO
from decorator import contextmanager
from sqlalchemy.exc import IntegrityError
from sqlalchemy.sql.functions import count

from sqlalchemy import bindparam, not_, literal
from sqlalchemy import and_
from sqlalchemy import func

from sqlalchemy.sql import select, exists
from sqlalchemy.sql import distinct
from sqlalchemy.sql.operators import in_op as in_

from pymlst.common import blat, psl, kma
from pymlst.common import utils, exceptions
from pymlst.wg import model


@contextmanager
def open_wg(file=None, ref='ref'):
    """A context manager function to wrap the creation a
       :class:`~pymlst.wg.core.WholeGenomeMLST` object.

    Context managers allow you to instantiate objects using the ``with`` keyword,
    eliminating the need to manage exceptions and commit/close processes yourself.

    :param file: The path to the database file to work with.
    :param ref: The name that will be given to the reference strain in the database.

    :yields: A :class:`~pymlst.wg.core.WholeGenomeMLST` object.
    """
    mlst = WholeGenomeMLST(file, ref)
    try:
        yield mlst
    finally:
        mlst.close()


class DuplicationHandling(Enum):
    CONCATENATE = 1
    REMOVE = 2


class DatabaseWG:
    """A core level class to manipulate the genomic database.

        .. warning:: Shouldn't be instantiated directly,
                     see :class:`~pymlst.wg.core.WholeGenomeMLST` instead.
    """

    def __init__(self, file, ref):
        """
        :param path: The path to the database file to work with.
        """
        engine = utils.get_updated_engine(file, 'wg')
        self.__connection = engine.connect()

        self.__cached_queries = {}

        self.__ref = ref
        self.__separator = ';'

        self.__core_genome = self.__load_core_genome()

    @contextmanager
    def begin(self):
        with self.__connection.begin():
            yield

    @property
    def ref(self):
        return self.__ref
    
    @property
    def separator(self):
        return self.__separator
    
    @property
    def connection(self):
        return self.__connection

    @property
    def core_genome(self):
        return self.__core_genome

    def __get_cached_query(self, name, query_supplier):
        if name in self.__cached_queries:
            return self.__cached_queries[name]
        query = query_supplier()
        self.__cached_queries[name] = query
        return query

    def __load_core_genome(self):
        result = self.connection.execute(
            select([model.mlst.c.gene, model.sequences.c.sequence])
            .where(and_(
                model.mlst.c.souche == self.__ref,
                model.mlst.c.seqid == model.sequences.c.id))
        ).fetchall()
        genes = {}
        for row in result:
            genes[row.gene] = row.sequence
        return genes

    def check_name(self, name, strain=False):
        if self.__separator in name:
            raise exceptions.InvalidGeneName(
                '{} contains {} symbol'.format(name, self.__separator))
        if strain:
            if "-" in name:
                logging.warning("Strain name '{}' contain '-', that could make some problems for further analysis".format(name))

    def add_core_genome(self, gene, sequence, mode=None):
        self.check_name(gene)
        if gene in self.__core_genome:
            raise exceptions.DuplicatedGeneName(
                '{} is duplicated'.format(gene))
        added, seq_id = self.__add_sequence(sequence)
        if not added:
            if mode is None:
                raise exceptions.DuplicatedGeneSequence(
                    'Duplicated sequence for gene {}'.format(gene))
            if mode == DuplicationHandling.CONCATENATE:
                self.__concatenate_gene(seq_id, gene)
            elif mode == DuplicationHandling.REMOVE:
                self.__remove_sequence(seq_id)
            return False
        self.__add_mlst(self.__ref, gene, seq_id)
        self.__core_genome[gene] = sequence
        return True

    def add_genome(self, gene, strain, sequence):
        _, seq_id = self.__add_sequence(sequence)
        self.__add_mlst(strain, gene, seq_id)

    def __add_mlst(self, souche, gene, seqid):
        """Adds an MLST gene bound to an existing sequence."""
        existing = self.connection.execute(
            select([model.mlst.c.id])
            .where(model.mlst.c.souche == souche, model.mlst.c.gene == gene, \
                   model.mlst.c.seqid == seqid)
        ).fetchone()
        if existing is None:
            self.connection.execute(
                model.mlst.insert()\
                .values(souche=souche, gene=gene, seqid=seqid))

    def __add_sequence(self, sequence):
        """Adds a sequence if it doesn't already exist."""
        query = self.__get_cached_query(
            'add_sequence',
            lambda:
            select([model.sequences.c.id])
                .where(model.sequences.c.sequence == bindparam('sequence')))

        existing = self.connection.execute(
            query,
            sequence=sequence
        ).fetchone()

        if existing is not None:
            return False, existing.id

        res = self.connection.execute(
            model.sequences.insert(),
            sequence=sequence)

        return True, res.inserted_primary_key[0]

    def __concatenate_gene(self, seq_id, gene_name):
        """Associates a new gene to an existing sequence using concatenation."""
        self.connection.execute(
            model.mlst.update()
                .values(gene=model.mlst.c.gene + ';' + gene_name)
                .where(model.mlst.c.seqid == seq_id))

    def __remove_sequence(self, seq_id):
        self.connection.execute(
            model.sequences.delete()
            .where(model.sequences.c.id == seq_id))
        self.connection.execute(
            model.mlst.delete()
            .where(model.mlst.c.seqid == seq_id))

    def __remove_orphan_sequences(self, ids):
        """Removes sequences if they aren't referenced by any gene."""
        query = model.sequences.delete() \
            .where(and_(
                not_(exists(
                    select([model.mlst.c.id])
                    .where(model.mlst.c.seqid == model.sequences.c.id))),
                model.sequences.c.id == bindparam('seqid')))

        for seqid in ids:
            self.connection.execute(
                query,
                seqid=seqid)

    def remove_gene(self, gene):
        """Removes a specific gene and its sequences."""
        ids = self.__get_gene_sequences_ids(gene)
        if len(ids) == 0:
            return False
        self.connection.execute(
            model.mlst.delete()
            .where(model.mlst.c.gene == gene))
        self.__remove_orphan_sequences(ids)
        if gene in self.core_genome:
            self.core_genome.pop(gene)
        return True

    def remove_strain(self, strain):
        """Removes a specific strain."""
        if strain == self.ref:
            raise exceptions.ReferenceStrainRemoval(
                '{} strain can not be removed'.format(self.ref))
        ids = self.__get_strain_sequences_ids(strain)
        if len(ids) == 0:
            return False
        self.connection.execute(
            model.mlst.delete()
                .where(model.mlst.c.souche == strain))
        self.__remove_orphan_sequences(ids)
        return True

    def __get_gene_sequences_ids(self, gene):
        """Return the IDs of the sequences associated with a specific gene."""
        rows = self.connection.execute(
            select([model.mlst.c.seqid])
            .where(model.mlst.c.gene == gene)
        ).fetchall()
        return {row.seqid for row in rows}

    def __get_strain_sequences_ids(self, strain):
        """Return the IDs of the sequences associated with a specific strain."""
        rows = self.connection.execute(
            select([model.mlst.c.seqid])
            .where(model.mlst.c.souche == strain)
        ).fetchall()
        return {row.seqid for row in rows}

    def contains_souche(self, souche):
        """Whether the strain exists in the base or not."""
        return self.connection.execute(
            select([literal(True)])
            .where(exists(
                select([model.mlst])
                .where(model.mlst.c.souche == souche)))
        ).scalar() is True  # -> True or False (Otherwise returns True or None)

    def get_gene_sequences(self, gene):
        """Return all the sequences for a specific gene and
        lists the strains that are referencing them."""
        query = self.__get_cached_query(
            'get_gene_sequences',
            lambda:
            select([model.mlst.c.seqid,
                    func.group_concat(model.mlst.c.souche, bindparam('separator')),
                    model.sequences.c.sequence])
            .select_from(
                model.mlst.join(model.sequences))
            .where(and_(
                model.mlst.c.gene == bindparam('gene'),
                model.mlst.c.souche != bindparam('souche')))
            .group_by(model.mlst.c.seqid))

        res = self.connection.execute(
            query,
            separator=';',
            souche=self.__ref,
            gene=gene
        ).fetchall()

        seqs = []
        for seq in res:
            tmp = seq[1].split(";")
            tmp.sort()
            seqs.append([seq[0], tmp, seq[2]])
        return seqs

    def get_gene_sequence_reference(self, gene):
        return(self.__core_genome.get(gene, []))
    
    def get_duplicated_genes(self):
        """Return the genes that are duplicated."""
        m_alias = model.mlst.alias()

        exist_sub = select([model.mlst]) \
            .where(and_(
                model.mlst.c.gene == m_alias.c.gene,
                model.mlst.c.souche == m_alias.c.souche,
                model.mlst.c.id != m_alias.c.id)) \
            .limit(1)

        res = self.connection.execute(
            select([model.mlst.c.gene])
            .where(and_(
                model.mlst.c.souche != self.__ref,
                exists(exist_sub)))
            .group_by(model.mlst.c.gene)
        ).fetchall()

        return {row[0] for row in res}

    def get_all_strains(self):
        """Return all distinct strains."""
        res = self.connection.execute(
            select([distinct(model.mlst.c.souche)]).
            where(model.mlst.c.souche != self.__ref)
        ).fetchall()
        return [r[0] for r in res]

    def get_core_genes(self):
        """Return all distinct genes."""
        res = self.connection.execute(
            select([distinct(model.mlst.c.gene)]).
            where(model.mlst.c.souche == self.__ref)
        ).fetchall()
        return [r[0] for r in res]

    def count_sequences_per_gene(self):
        """Return the number of distinct sequences per gene."""
        res = self.connection.execute(
            select([model.mlst.c.gene, count(distinct(model.mlst.c.seqid))])
            .where(model.mlst.c.souche != self.__ref)
            .group_by(model.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_souches_per_gene(self):
        """Return the number of distinct stains per gene."""
        res = self.connection.execute(
            select([model.mlst.c.gene, count(distinct(model.mlst.c.souche))])
            .where(model.mlst.c.souche != self.__ref)
            .group_by(model.mlst.c.gene)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_genes_per_souche(self, valid_shema):
        """Return the number of distinct genes per strain.

        The counted genes are restricted to the ones given in the valid_schema."""
        res = self.connection.execute(
            select([model.mlst.c.souche, count(distinct(model.mlst.c.gene))])
            .where(in_(model.mlst.c.gene, valid_shema))
            .group_by(model.mlst.c.souche)
        ).fetchall()
        return {r[0]: r[1] for r in res}

    def count_sequences(self):
        """Return the number of distinct."""
        return self.connection.execute(
               select([count(distinct(model.mlst.c.seqid))])
               .where(model.mlst.c.souche != self.__ref)
        ).fetchone()[0]

    def get_strains_distances(self, valid_schema):
        """Return the distances between strains.

        For all the possible pairs of strains, counts how many of their genes
        are different (different seqids so different sequences).
        The compared genes are restricted to the ones given in the valid_schema.
        """
        # alias_1 = model.mlst.alias()
        # alias_2 = model.mlst.alias()

        # result = self.connection.execute(
        #     select(
        #         [alias_1.c.souche, alias_2.c.souche, count(distinct(alias_1.c.gene))])
        #     .select_from(
        #         alias_1.join(
        #             alias_2,
        #             and_(
        #                 alias_1.c.seqid != alias_2.c.seqid,
        #                 alias_1.c.souche != alias_2.c.souche,
        #                 alias_1.c.gene == alias_2.c.gene)))
        #     .where(
        #         and_(
        #             in_(alias_1.c.gene, valid_schema),
        #             alias_1.c.souche != self.__ref,
        #             alias_2.c.souche != self.__ref))
        #     .group_by(alias_1.c.souche, alias_2.c.souche)
        # ).fetchall()

        # distance = {}
        # for entry in result:
        #     dist = distance.setdefault(entry[0], {})
        #     dist[entry[1]] = entry[2]

        mlst = self.get_mlst(valid_schema)
        strains = self.get_all_strains()
        
        logging.info("Start distance calculation. This can take while")
        distance = {}
        
        for i,s1 in enumerate(strains):
            dist = distance.setdefault(s1, {})
            for s2 in strains:
                count = 0
                if s1==s2:
                    dist[s2] = 0
                    continue
                for k in mlst.values():
                    if s1 in k and s2 in k:
                        if k.get(s1) != k.get(s2):
                            count +=1
                dist[s2] = count
            logging.info("Strains %s on "+str(len(strains)), str(i))
        
        return distance

    def get_mlst(self, valid_schema):
        """Return the the genes MLST.

        Returns a dictionary associated with each gene, all the strains that
        reference it, and their sequence ids. The genes returned are restricted
        to those given in the valid_schema.
        """
        result = self.connection.execute(
            select([model.mlst.c.gene, model.mlst.c.souche,
                    func.group_concat(model.mlst.c.seqid, ';')])
            .where(and_(model.mlst.c.souche != self.__ref,
                        in_(model.mlst.c.gene, valid_schema)))
            .group_by(model.mlst.c.gene, model.mlst.c.souche)
        ).fetchall()

        mlst = {}

        for entry in result:
            sequences = mlst.setdefault(entry[0], {})
            sequences[entry[1]] = entry[2]
        return mlst

    def close(self):
        self.__connection.close()


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
        self.__database = DatabaseWG(file, ref)
        self.__file = file
        utils.create_logger()

    @property
    def database(self):
        return self.__database

    def create(self, coregene, concatenate=False, remove=False):
        """Creates a whole genome MLST database from a core genome `fasta`_ file.

        :param coregene: The fasta file containing the reference core genome.
        :param concatenate: Whether we should concatenate genes with identical sequences.
        :param remove: Whether we should remove genes with identical sequences.

        For instance, if concatenate is set to ``True``, 2 genes **g1** and **g2**
        having the same sequence
        will be stored as a single gene named **g1;g2**.
        """
        ##remove old indexing
        kma.delete_indexing(self.__file)
        
        with self.database.begin():
            rc_genes = 0
            invalid_genes = 0

            for gene in SeqIO.parse(coregene, 'fasta'):
                if not utils.validate_sequence(gene.seq):
                    gene.seq = gene.seq.reverse_complement()
                    if utils.validate_sequence(gene.seq):
                        rc_genes += 1
                    else:
                        invalid_genes += 1
                        logging.debug("Skipped Invalid gene %s", gene.id)
                        continue

                if concatenate:
                    mode = DuplicationHandling.CONCATENATE
                elif remove:
                    mode = DuplicationHandling.REMOVE
                else:
                    mode = None

                ##clean geneid
                geneid = utils.clean_geneid(gene.id)
                added = self.__database.add_core_genome(geneid, str(gene.seq), mode)

                if not added:
                    if concatenate:
                        logging.debug("Concatenated gene %s", gene.id)
                    elif remove:
                        logging.debug('Gene %s has duplicated sequence, removed', gene.id)

            if rc_genes:
                logging.info('Reverse-complemented genes: %s', str(rc_genes))

            if invalid_genes:
                logging.info('Skipped invalid genes: %s', str(invalid_genes))

            if len(self.__database.get_core_genes()) == 0:
                raise exceptions.InvalidGeneName("No valid gene found\n" + \
                    "You probably load genome instead of genes")

            logging.info('Database initialized')

    def add_strain(self, genome, strain=None, identity=0.95, coverage=0.90):
        """Adds a genome strain to the database.

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
        with self.database.begin():
            if identity < 0 or identity > 1:
                raise exceptions.BadIdentityRange('Identity must be in range [0-1]')
            if coverage <0 or coverage > 1:
                raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')

            name = strain
            if name is None:
                name = genome.name.split('/')[-1]
            self.__database.check_name(name, strain=True)

            tmpfile, tmpout = blat.blat_tmp()
            tmpout.close()

            try:
                # verify that the strain is not already in the database
                if self.__database.contains_souche(name):
                    raise exceptions.StrainAlreadyPresent(
                        'Strain {} already present in the base'.format(name))

                # read coregene
                self.__create_core_genome_file(tmpfile)
                tmpfile.close()

                # BLAT analysis
                logging.info("Search coregene with BLAT")
                genes = blat.run_blat(genome, tmpfile, tmpout, identity, coverage)
                logging.info("Finish run BLAT, found %s genes", len(genes))

                # add MLST sequence
                seqs = utils.read_genome(genome)

                bad = 0
                partial = 0
                partial_filled = 0

                for core_gene in self.__database.core_genome.keys():

                    if core_gene not in genes:
                        continue

                    for gene in genes.get(core_gene):

                        seq = seqs.get(gene.chro, None)
                        if seq is None:
                            raise exceptions.ChromosomeNotFound(
                                'Chromosome {} not found '.format(gene.chro))

                        # Correct coverage
                        if gene.coverage == 1:
                            sequence = gene.get_sequence(seq)
                        else:
                            coregene_seq = self.__database.core_genome[core_gene]
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
                        self.__database.add_genome(gene.gene_id(), name, str(sequence))

                logging.info("Added %s new MLST genes to the database", len(genes) - bad)
                logging.info('Found %s partial genes, filled %s', partial, partial_filled)
                logging.info('Removed %s genes', bad)
                logging.info("DONE")

            finally:
                if os.path.exists(tmpfile.name):
                    os.remove(tmpfile.name)
                if os.path.exists(tmpout.name):
                    os.remove(tmpout.name)

    def add_reads(self, fastqs, strain=None, identity=0.95, coverage=0.90, \
                  reads=10):
        """
        Adds raw reads of a strain to the database.

        How it works:

        1. A `KMA`_ research is performed on reads (fastq) of the strain to find
           sub-sequences matching the core genes.
        2. The identified sub-sequences are extracted and added to our database where
           they are associated to a **sequence ID**.
        3. An MLST entry is created, referencing the sequence,
           the gene it belongs to, and the strain it was found in.

        :param fastqs: The reads we want to add as a list of `fastq`_ file.
        :param strain: The name that will be given to the new strain in the database.
        :param identity: Sets the minimum identity used by `BWA`_ for sequences research (in percent).
        :param coverage: Sets the minimum accepted coverage for found sequences.
        :param reads: Sets the minimum number of reads coverage to conserved an results 
        """
        with self.database.begin():
            if identity < 0 or identity > 1:
                raise exceptions.BadIdentityRange('Identity must be in range [0-1]')
            if coverage <0 or coverage > 1:
                raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')

            ##indexing
            if kma.is_database_indexing(self.__file) is False:
                with kma.index_tmpfile() as tmpfile:
                    coregene = self.__create_core_genome_file(tmpfile)
                    tmpfile.flush()
                    kma.index_database(self.__file, tmpfile)

            ##Strain name
            name = strain
            if name is None:
                name = fastqs[0].name.split('/')[-1]
            self.__database.check_name(name)
        
            ##run kma
            kma_res,seqs = kma.run_kma(fastqs, self.__file, identity, coverage, reads)
            core_genes = self.__database.core_genome
            
            valid = 0
            minus = 0
            frame = 0
            for res in kma_res:
                seq = seqs.get(res)
                if seq is None:
                    raise PyMLSTError("%s not found in the fasta files", res)

                ## test minus
                b = (seq.count('a') + seq.count('t') + seq.count('c') + \
                     seq.count('g'))
                if b != 0:
                    minus +=1
                    logging.debug("%s Remove incertain", res)
                    continue

                ## test CDS
                try:
                    seq.translate(cds=True, table=11)
                except:
                    frame += 1
                    logging.debug("%s Remove bad CDS", res)
                    continue
 
                ##add sequence and MLST
                gene = res.split("_")[0]
                if gene not in core_genes:
                    logging.warning("Gene %s not present in database", gene)
                    continue
                valid +=1
                self.__database.add_genome(gene, name, str(seq))
                
            logging.info("Add %s new MLST genes to database", str(valid))
            logging.info("Remove %s genes with uncertain bases", str(minus))
            logging.info("Remove %s genes with bad CDS", str(frame))  

            
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
            raise exceptions.NothingToRemove('No gene to remove found')
        all_genes = set(all_genes)

        for gene in all_genes:
            if self.__database.remove_gene(gene):
                logging.info("%s : OK", gene)
            else:
                logging.warning("%s : Not found", gene)

        ##delete kma indexing as gene have been modified
        if kma.is_database_indexing(self.__file):
            kma.delete_indexing(self.__file)

    def remove_strain(self, strains, file=None):
        """Removes entire strains from the database.

        :param strains: Names of the strains to remove.
        :param file: A file containing a strain name per line.
        """
        if self.__database.ref in strains:
            raise exceptions.ReferenceStrainRemoval(
                '{} strain can not be removed'.format(self.__database.ref))

        # list strains to remove
        all_strains = utils.strip_file(file)
        if strains is not None:
            all_strains.extend(strains)
        if len(all_strains) == 0:
            raise exceptions.NothingToRemove('No strain to remove found')
        all_strains = set(all_strains)

        for strain in all_strains:
            if self.__database.remove_strain(strain):
                logging.info("%s : OK", strain)
            else:
                logging.warning("%s : Not found", strain)

    def extract(self, extractor, output=sys.stdout):
        """Takes an extractor object and writes the extraction result on the given output.

        :param extractor: A :class:`~pymlst.wg.core.Extractor` object describing
                          the way data should be extracted.
        :param output: The output that will receive extracted data.
        """
        extractor.extract(self.__database, output)

    def __create_core_genome_file(self, tmp_file):
        ref_genes = self.__database.core_genome
        for gene, sequence in ref_genes.items():
            tmp_file.write('>' + gene + "\n" + sequence + "\n")

    def close(self):
        self.database.close()


class Extractor(ABC):
    """A simple interface to ease the process of creating new extractors."""
    @abstractmethod
    def extract(self, base, output):
        """
        :param base: The database to extract data from.
        :param output: The output where to write the extraction results.
        """


def find_recombination(genes, alignment, output=sys.stdout):
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
            raise exceptions.PyMLSTError(
                'The genes list doesn\'t correspond to the alignment {}'.format(indice))

        # genes
        sequences[indice].append(line)
        indice += 1

    # check sequences are correctly align
    for i, seqs in enumerate(sequences):
        if len({len(s) for s in seqs}) > 1:
            logging.error({len(s) for s in seqs})
            raise exceptions.PyMLSTError(
                'The following genes are not aligned: {}'.format(genes[i]))

    output.write("Gene\tMutation\tLenght\tmutation per 100 base\n")
    for i, seqs in enumerate(sequences):
        compared = utils.compar_seqs(seqs)
        output.write(genes[i] + "\t" + str(compared) + "\t" + str(len(seqs[0])) + \
                     "\t" + str(compared/len(seqs[0])*100) + "\n")


def find_subgraph(distance, threshold=50, output=sys.stdout, export='list'):
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
        raise exceptions.PyMLSTError(
            "The distance file seems not correctly "
            "formatted, not integer on first line") from err

    for line in distance.readlines():
        dist_line = line.rstrip("\n").split("\t")
        samps.append(dist_line[0])
        dists.append(dist_line[1:])

    if len(samps) != strains:
        raise exceptions.PyMLSTError(
            "The distance is not properly formatted, "
            "the number of strains ({}) doesn't correspond to {}"
            .format(len(samps), strains))

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
