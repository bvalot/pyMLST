"""Core classes and functions to work with Classical MLST data."""

from io import TextIOWrapper
import logging
import os
import re
import sys

from Bio import SeqIO
from alembic.command import heads
from decorator import contextmanager

from sqlalchemy import and_
from sqlalchemy.exc import IntegrityError

from sqlalchemy.sql import select
from sqlalchemy.sql import distinct

from pymlst.cla import model
from pymlst.common import blat, kma,  utils, exceptions
from pymlst.common.utils import read_genome, create_logger


@contextmanager
def open_cla(file=None, ref=1):
    """A context manager function to wrap the creation a
       :class:`~pymlst.cla.core.ClassicalMLST` object.

    Context managers allow you to instantiate objects using the ``with`` keyword,
    eliminating the need to manage exceptions and commit/close processes yourself.

    :param file: The path to the database file to work with.
    :param ref: The name that will be given to the reference strain in the database.

    :yields: A :class:`~pymlst.cla.core.ClassicalMLST` object.
    """
    mlst = ClassicalMLST(file, ref)
    try:
        yield mlst
    finally:
        mlst.close()


class DatabaseCLA:
    """A core level class to manipulate the genomic database.

        .. warning:: Shouldn't be instantiated directly,
                     see :class:`~pymlst.cla.core.ClassicalMLST` instead.
    """

    def __init__(self, file, ref):
        """
        :param path: The path to the database file to work with.
        """
        engine = utils.get_updated_engine(file, 'cla')
        self.__connection = engine.connect()

        self.__ref = ref

        self.__core_genome = self.__load_core_genome()

    @contextmanager
    def begin(self):
        with self.__connection.begin():
            yield

    @property
    def ref(self):
        return self.__ref

    @property
    def connection(self):
        return self.__connection

    @property
    def core_genome(self):
        return self.__core_genome

    def __load_core_genome(self):
        return self.get_genes_by_allele(self.__ref)

    def add_sequence(self, sequence, gene, allele):
        """Adds a new sequence associated to a gene and an allele."""
        try:
            self.connection.execute(
                model.sequences.insert(),
                sequence=sequence, gene=gene, allele=allele)
        except IntegrityError:
            raise exceptions.DuplicatedGeneSequence(
                'Duplicated sequence for gene {} and allele {}'.format(gene, allele))
        if allele == self.__ref:
            if gene not in self.__core_genome:
                self.__core_genome[gene] = sequence

    def add_mlst(self, sequence_typing, gene, allele):
        """Adds a new sequence typing, associated to a gene and an allele."""
        # sequence = self.get_sequence_by_gene_and_allele(gene, allele)
        # if sequence is None:
        #     raise exceptions.AlleleSequenceNotFound(
        #         'No sequence found for allele {} of gene {}'
        #         .format(allele, gene))
        self.connection.execute(
            model.mlst.insert(),
            st=sequence_typing, gene=gene, allele=allele)

    def get_genes_by_allele(self, allele):
        """Returns all the distinct genes in the database and their sequences for a given allele."""
        genes = self.connection.execute(
            select([distinct(model.mlst.c.gene), model.sequences.c.sequence])
            .select_from(model.mlst.join(
                model.sequences,
                model.mlst.c.gene == model.sequences.c.gene))
            .where(model.sequences.c.allele == allele)
        ).fetchall()
        return {gene.gene: gene.sequence for gene in genes}

    def get_allele_by_sequence_and_gene(self, sequence, gene):
        """Gets an allele by sequence and gene."""
        res = self.connection.execute(
            select([model.sequences.c.allele])
            .where(and_(
                model.sequences.c.sequence == sequence,
                model.sequences.c.gene == gene))
        ).fetchone()
        if res is None:
            return None
        return res.allele

    def get_st_by_gene_and_allele(self, gene, allele):
        """Gets all the STs of a gene/allele pair."""
        seq_types = self.connection.execute(
            select([model.mlst.c.st])
            .where(and_(
                model.mlst.c.gene == gene,
                model.mlst.c.allele == allele))
        ).fetchall()
        return [seq_t.st for seq_t in seq_types]

    def get_sequence_by_gene_and_allele(self, gene, allele):
        """Gets a sequence by gene and allele."""
        res = self.connection.execute(
            select([model.sequences.c.sequence])
            .where(and_(
                model.sequences.c.gene == gene,
                model.sequences.c.allele == allele))
        ).fetchone()
        if res is None:
            return None
        return res.sequence

    def get_all_sequences_by_gene(self, gene):
        """Get all the sequences from a particular gene"""
        seqs = self.connection.execute(
            select([model.sequences.c.allele, model.sequences.c.sequence])
            .where(model.sequences.c.gene == gene)
        ).fetchall()
        return {seq.sequence:seq.allele for seq in seqs}     

    def get_all_sequences(self):
        """Get all sequences allele"""
        seqs = self.connection.execute(
            select([model.sequences.c.allele, model.sequences.c.gene, model.sequences.c.sequence])
        ).fetchall()
        return {seq.gene+"_"+str(seq.allele):seq.sequence for seq in seqs}

    def remove_allele(self, gene, allele):
        """Remove an particular allele for the gene"""
        res1 = self.connection.execute(
            select([model.sequences.c.id])
            .where(and_(model.sequences.c.gene == gene,
                        model.sequences.c.allele == allele))
            ).fetchone()
        if res1 is None:
            logging.warning("No sequence found for the gene %s with the allele %s", gene, str(allele))
            return(False)
        else:
            self.connection.execute(
                model.sequences.delete()
                .where(model.sequences.c.id == res1.id))
            logging.debug("Delete sequence id %s", str(res1.id))
                
        ##Get ST number with this allele and remove its
        res2 =  self.connection.execute(
            select([model.mlst.c.st])
            .where(and_(
                model.mlst.c.gene == gene,
                model.mlst.c.allele == allele))
        ).fetchall()
        for res in res2:
            continue
            self.connection.execute(
                model.mlst.delete()
                .where(model.mlst.c.st == res.st))
        logging.info("Remove %s ST containing this allele", str(len(res2)))
        return(True)
    
    def close(self):
        self.__connection.close()


class ClassicalMLST:
    """Classical MLST python representation.

        Example of usage::

            open_cla('database.db') as db:
                db.create(open('scheme.txt'), [open('gene1.fasta'),
                                               open('gene2.fasta'),
                                               open('gene3.fasta')])
                db.multi_search(open('genome.fasta'))
        """

    def __init__(self, file, ref):
        """
        :param file: The path to the database file to work with.
        :param ref: The name that will be given to the reference strain in the database.
        """
        self.__database = DatabaseCLA(file, ref)
        self.__file = file
        create_logger()

    @property
    def database(self):
        return self.__database

    def create(self, scheme, alleles):
        """Creates a classical MLST database from an MLST profile and a list of alleles.

        :param scheme: The MLST profile
        :param alleles: A list of alleles files.

        The MLST profile should be a **CSV** file respecting the following format:

        .. csv-table:: MLST Profile CSV
            :header: "ST", "gene1", "gene2", "gene3", "..."
            :widths: 5, 5, 5, 5, 5

            1, 1, 1, 1, ...
            2, 3, 3, 2, ...
            3, 1, 2, 1, ...
            4, 1, 1, 3, ...
            ..., ..., ..., ..., ...

        """
        ##remove old indexing
        kma.delete_indexing(self.__file)
        
        with self.database.begin():
            # Verify sheme list with fasta files
            header = scheme.readline().rstrip("\n").split("\t")
            if len(header) != len(alleles) + 1:
                raise Exception("The number of genes in sheme don't "
                                "correspond to the number of fasta file\n"
                                + " ".join(header) + "\n")
            fastas = {}
            for fi in alleles:
                name = fi.name.split("/")[-1]
                name = name[:name.rfind(".")]
                if name not in header:
                    raise Exception("Gene " + name + " not found in sheme\n" + " ".join(header))
                fastas[name] = fi

            # load sequence allele
            alleles = {}
            for gene, fi in fastas.items():
                logging.debug("Read sequence from %s", fi.name)
                alleles[gene] = set()
                for seq in SeqIO.parse(fi, 'fasta'):
                    try:
                        match = re.search('[0-9]+$', seq.id)
                        allele = int(match.group(0))
                    except Exception as err:
                        raise Exception("Unable to obtain allele number for the sequence: " + seq.id) from err

                    self.__database.add_sequence(str(seq.seq), gene, allele)
                    alleles.get(gene).add(allele)

            # load MLST sheme
            logging.debug("Load schema %s", scheme.name)
            for line in scheme:
                line_content = line.rstrip("\n").split("\t")
                sequence_typing = int(line_content[0])
                for gene, allele in zip(header[1:], line_content[1:]):
                    if int(allele) not in alleles.get(gene):
                        logging.warning(
                            "Unable to find the allele number %s"
                            " for gene %s; skipped", str(allele), gene)
                        ##self.__database.add_mlst(sequence_typing, gene, 0)
                    else:
                        self.__database.add_mlst(sequence_typing, gene, int(allele))

            logging.info('Database initialized')

    def search_st(self, genome, identity=0.90, coverage=0.90, fasta=None):
        """Search the **Sequence Type** number of a strain.

        :param genome: The strain genome we want to add as a `fasta`_ file.
        :param identity: Sets the minimum identity used by `BLAT`_
                         for sequences research (in percent).
        :param coverage: Sets the minimum accepted gene coverage for found sequences.
        :param fasta: A file where to export genes alleles results in a fasta format.
        """

        if identity < 0 or identity > 1:
            raise exceptions.BadIdentityRange("Identity must be between 0 to 1")
        if coverage <0 or coverage > 1:
            raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')

        tmpfile, tmpout = blat.blat_tmp()

        try:
            # read coregene
            self.__create_core_genome_file(tmpfile)
            core_genome = self.__database.core_genome
            tmpfile.close()

            sorted_genes = [gene for gene in core_genome]
            sorted_genes.sort()
            genome_name = os.path.basename(genome.name).split('.')[0]

            # BLAT analysis
            #logging.info("Search coregene with BLAT")
            allele = {i: [] for i in core_genome}
            logging.info("Search coregene with BLAT for %s", str(genome_name))
            try:
                genes = blat.run_blat(genome, tmpfile, tmpout, identity, coverage)
            except exceptions.CoreGenomePathNotFound:
                logging.warning("No gene found for the strain %s", genome_name)
                return ST_result(genome_name, [], allele)
            logging.info("Finish run BLAT, found %s genes", str(len(genes)))

            # Search sequence MLST
            seqs = read_genome(genome)
            logging.info("Search allele gene to database")
            # print(genes)
            sequence_type = {i: set() for i in core_genome}
            for core_gene, core_seq in core_genome.items():
                if core_gene not in genes:
                    allele.get(core_gene).append("")
                    continue
                for gene in genes.get(core_gene):
                    seq = seqs.get(gene.chro, None)
                    if seq is None:
                        raise Exception("Chromosome ID not found " + gene.chro)

                    # verify coverage and correct
                    # if gene.coverage != 1:
                    #     gene.searchCorrect()
                    #     logging.info("Gene %s fill: added", gene.gene_id())

                    if gene.coverage == 1:
                        sequence = gene.get_sequence(seq)
                    else:
                        logging.debug("Align incomplet gene : %s", core_gene)
                        sequence = gene.get_aligned_sequence(seq, core_seq)

                    sequence = str(sequence)

                    # write fasta file with coregene
                    if fasta is not None:
                        fasta.write(">" + genome_name + "|" + core_gene + "\n" )
                        fasta.write(sequence + "\n")

                    # search allele
                    res = self.__database.get_allele_by_sequence_and_gene(sequence, core_gene)
                    if res is not None:
                        allele.get(core_gene).append(str(res))
                        # cursor.execute('''SELECT st FROM mlst WHERE gene=? and allele=?''',
                        #                (coregene, row[0]))
                        seq_types = self.__database.get_st_by_gene_and_allele(core_gene, res)
                        for seq_type in seq_types:
                            sequence_type.get(core_gene).add(seq_type)
                    else:
                        allele.get(gene.gene_id()).append("new")

            # if only know allele or not found
            # Seach st
            # st_val = []
            # if sum([len(i) == 1 and i[0] != "new" for i in allele.values()]) == len(allele):
            #     tmp = None
            #     for st_value in sequence_type.values():
            #         if st_value:
            #             if tmp is None:
            #                 tmp = st_value
            #             else:
            #                 tmp = tmp.intersection(st_value)
            #     st_val = list(tmp)
            st_val = self.__determined_ST(allele, sequence_type)
                
            return ST_result(genome_name, st_val, allele)
        
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
            if os.path.exists(tmpout.name):
                os.remove(tmpout.name)
                
                
    def multi_search(self, genomes, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        """Searches the **Sequence Type** number of one or multi strain(s) from an assembly genome.
            
        :param genomes: Tuple of one or more strain genomes given as input
        :param output: An output for the sequence type research results.
        :param identity: Sets the minimum identity used by `BLAT`_
                         for sequences research (in percent).
        :param fasta: A file where to export genes alleles results in a fasta format.
        :param coverage: Sets the minimum accepted coverage for found sequences.
        """
        header = True
        for genome in genomes:
            logging.info("Search ST from files: %s", os.path.basename(genome.name)) 
            res = self.search_st(genome, identity, coverage, fasta)
            res.write(output, header)
            if header:
                header=False
            logging.info("FINISH")

    def search_read(self, fastqs, identity=0.9, coverage=0.95, reads=10, fasta=None):
        """Searches the **Sequence Type** from raw reads of one strain.

        :param fastq: List of fastq files containing raw reads
        :param identity: Sets the minimum identity used by `KMA`_
                         for sequences research (in percent).
        :param coverage: Sets the minimum accepted gene coverage for found sequences.
        :param reads: Sets the minimum reads coverage to conserve an mapping 
        :param fasta: A file where to export genes alleles results in a fasta format.
        """
        if identity < 0 or identity > 1:
            raise exceptions.BadIdentityRange("Identity must be between 0 to 1")
        if coverage <0 or coverage > 1:
            raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')
        
        ##indexing
        if kma.is_database_indexing(self.__file) is False:
            with kma.index_tmpfile() as tmpfile:
                self.__create_all_core_genome_file(tmpfile)
                tmpfile.flush()
                kma.index_database(self.__file, tmpfile)
                
        ##run kma
        genome_name = os.path.basename(fastqs[0].name).split('.')[0]
        core_genome = self.__database.core_genome
        allele = {i: [] for i in core_genome}
        try:
            kma_res,seqs = kma.run_kma(fastqs, self.__file, identity, coverage, reads)
        except exceptions.CoreGenomePathNotFound:
            logging.warning("No gene found for the strain %s", genome_name)
            return ST_result(genome_name, [], allele)
        
        logging.info("Search allele gene to database")
        sequence_type = {i: set() for i in core_genome}
        for res in kma_res:
            logging.debug("Looking for kma result of gene %s", res)
            sequence = seqs.get(res)
            if sequence is None:
                raise PyMLSTError("%s not found in the fasta files", res)

            ## test minus
            b = (sequence.count('a') + sequence.count('t') + sequence.count('c') + \
                 sequence.count('g'))
            if b != 0:
                logging.debug("%s Remove incertain", res)
                continue
            
            ##add sequence and MLST
            gene = res.split("_")[0]
            if gene not in core_genome:
                logging.warnings("Gene %s not present in database", gene)
                continue
            alle = self.__database.get_allele_by_sequence_and_gene(str(sequence), gene)
        
            # write fasta file with coregene
            if fasta is not None:
                fasta.write(">" + genome_name + "|" + gene + "\n" )
                fasta.write(str(sequence) + "\n")

            if alle is not None:
                allele.get(gene).append(str(alle))
                seq_types = self.__database.get_st_by_gene_and_allele(gene, alle)
                for seq_type in seq_types:
                    sequence_type.get(gene).add(seq_type)
            else:
                ##search substring sequence
                sequence_gene = self.__database.get_all_sequences_by_gene(gene)
                find = False
                for s in sequence_gene.keys() :
                    if s.find(str(sequence))!=-1:
                        find = True
                        alle = sequence_gene.get(s)
                        logging.debug("Find subsequence of gene %s with the allele %s", gene, str(alle))
                        allele.get(gene).append(str(alle))
                        seq_types = self.__database.get_st_by_gene_and_allele(gene, alle)
                        for seq_type in seq_types:
                            sequence_type.get(gene).add(seq_type)
                if find is False:
                    allele.get(gene).append("new")
        
        st_val = self.__determined_ST(allele, sequence_type)
        return ST_result(genome_name, st_val, allele)   

    def multi_read(self, fastqs, identity=0.90, coverage=0.95, reads=10, paired=True, fasta=None, output=sys.stdout):
        
        """Search the **Sequence Type** number of one or multi strain(s) from raw reads.
            
        :param fastqs: Tuple of one or more strain raw reads given as input
        :param output: An output for the sequence type research results.
        :param identity: Sets the minimum identity used by `KMA`_
                         for sequences research (in percent).
        :param reads: Sets the minimum reads coverage to conserve an mapping
        :param paired: Defined if the raxw reads are by paired or single
        :param fasta: A file where to export genes alleles results in a fasta format.
        :param coverage: Sets the minimum accepted coverage for found sequences.
        
        """
        header = True
        if paired:
            if len(fastqs)%2 != 0:
                raise exceptions.PyMLSTError("Fastq paired files are not a multiple of 2")
            for fastq in zip(fastqs[::2],fastqs[1::2]):
                logging.info("Search ST from files: %s - %s", os.path.basename(fastq[0].name), \
                             os.path.basename(fastq[1].name))
                res = self.search_read(fastq, identity, coverage, reads, fasta)
                res.write(output, header)
                if header:
                    header=False
                logging.info("FINISH")
        else:
            for fastq in fastqs:
                logging.info("Search ST from files: %s", os.path.basename(fastq.name)) 
                res = self.search_read([fastq], identity, coverage, reads, fasta)
                res.write(output, header)
                if header:
                    header=False
                logging.info("FINISH")

    def remove_allele(self, gene, allele):
        """Removes some problematic allele out of the claMLST database

        :param gene: Gene name on the database, ignoring case
        :param allele: Integer allele number of this gene
        """
        if self.__database.remove_allele(gene, allele):
            kma.delete_indexing(self.__file)
        

    def __create_core_genome_file(self, tmpfile):
        core_genome = self.__database.core_genome
        for gene, sequence in core_genome.items():
            tmpfile.write('>' + gene + "\n" + sequence + "\n")

    def __create_all_core_genome_file(self, tmpfile):
        for gene, sequence in self.__database.get_all_sequences().items():
            tmpfile.write('>' + gene + "\n" + sequence + "\n")

    def __determined_ST(self, allele, sequence_type):
        # Seach st
        st_val = []
        if sum([len(i) == 1 and i[0] != "new" for i in allele.values()]) == len(allele):
            tmp = None
            for st_value in sequence_type.values():
                if st_value:
                    if tmp is None:
                        tmp = st_value
                    else:
                        tmp = tmp.intersection(st_value)
            st_val = list(tmp)
        return st_val
        
    def close(self):
        self.__database.close()


class ST_result:
    """ Writes the results of the ST research"""
    
    def __init__(self,genome_name, st_val, alleles):
        """
        :param genome_name: Name of the genome retrieved from the path provided by the user
        :param st_val:  ST values identified for each genome by search_st
        :param alleles: Alleles of the strain recovered in the core genome
        """
        self.name = genome_name
        self.st = st_val
        self.alleles = alleles
        
    def write(self,output=sys.stdout, header=True):
        """Writes the results in output file"""
        genes = self.sort_genes()
        if header:
            output.write("Sample\tST\t" + "\t".join(genes) + "\n")
        towrite = [self.name]
        towrite.append(self.__str_list(self.st))
        towrite.extend([self.__str_list(self.alleles.get(g, "")) for g in genes])
        output.write("\t".join(towrite) + "\n")
    
    def sort_genes(self):
        genes =  list(self.alleles.keys())
        genes.sort()
        return genes
        
    def __str_list(self,l):
        return ";".join(map(str, l))
