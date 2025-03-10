""" Core classes and functions to work in alternative typing methods. """

import logging
import os
import sys
import re
import io
import tempfile


from alembic.command import heads
from Bio import SeqIO
from decorator import contextmanager


from sqlalchemy import and_
from sqlalchemy.exc import IntegrityError

from sqlalchemy.sql import select
from sqlalchemy.sql import distinct

# For abstract class PyTyper
from abc import ABC, abstractmethod

from pymlst.pytyper import model
from pymlst.pytyper.method import FIM, SPA, CLMT
from pymlst.pytyper.url import SPA_URL_TYPE, SPA_URL_SEQ, FIM_URL
from pymlst import config
from pymlst.common import blat, kma, utils, exceptions, web
from pymlst.common.utils import create_logger, read_genome, file_name

# Creates generator function for pytyper objects (and add context manager decorator) 
@contextmanager
def open_typer(method):
    """
    :param method: Defines typing method to apply. Possible values :
        1- fim
        2- spa
        3- clmt

    :yields: A :class: '~pymlst.pytyper.core.pyTyper' object.
    """
    utils.create_logger()
    
    # absolute path to automatically initialized pyTyper database
    fi = config.get_data('pytyper.db')

    if method == FIM:
        typer = FimH(fi)
    elif method == SPA:
        typer = Spa(fi)
    elif method == CLMT:
        typer = Clmt(fi)
    else:
        raise exceptions.WrongBaseType(f'Method {method} not implemented.') 
    
    try:
        yield typer
    finally:
        typer.close()



class DatabaseTyper:

    def __init__(self, fi):
        """
        :param file: The path to the database file to work with.
        """
        engine = utils.get_updated_engine(fi, 'pytyper')
        self.__connection = engine.connect()
        self.__file = fi

    @contextmanager
    def begin(self):
        with self.__connection.begin():
            yield

    @property
    def connection(self):
        return self.__connection

    def check_db(self, method):
        """Checks if database contains specified typing method."""
        count = self.connection.execute(
                select([model.typerSeq.c.id])
                .where(model.typerSeq.c.typing == method)
                ).fetchall()

        return(len(count) > 1)


    def add_sequence(self, sequence, method, allele):
        try:
            self.connection.execute(
                model.typerSeq.insert(),
                sequence=sequence, typing=method, allele=allele
                    )
        except IntegrityError:
            logging.warning(f"Duplicate sequence for allele {allele}.")
            return(False)
        return(True)

    def get_sequences(self, method):
        """Gets all sequences for an method"""
        seqs = self.connection.execute(
            select([model.typerSeq.c.allele, model.typerSeq.c.sequence])
            .where(model.typerSeq.c.typing == method)
        ).fetchall()
        return(seqs)
        
    def add_st(self, st, method, allele):
        self.connection.execute(
            model.typerSt.insert(),
            st=st, typing=method, allele=allele
        )

    def get_st(self, method, allele):
        """Gets the STs of an allele."""
        typer_st = self.connection.execute(
            select(model.typerSt.c.st)
            .where(and_(
                model.typerSt.c.typing == method,
                model.typerSt.c.allele == allele))
        ).first()
        if typer_st is None or typer_st[0] == 'NT':
            return ""
        else:
            return(typer_st[0])
        
    def close(self) -> None:
        self.__connection.close()


# Generic class for pyTyper
class PyTyper(ABC):
    """ Primary class for all pyTyper objects listed on method
    """
    def __init__(self, fi , typing):
        """ 
        :param fi: Path to the database file
        :param typing: Typing type
        """
        self._database = DatabaseTyper(fi)
        self._typing = typing
        if not self._database.check_db(self._typing):
            self.create()
        

    @abstractmethod
    def search_genome(self, genome, identity=0.90, coverage=0.90, fasta=None):
        """
        Abstract method for searching alleles against a genome.

        :param genome: Path to the fasta genome
        :param identity: Minimum identity treshold (0.9)
        :param coverage: Minimum coverage threshold (0.9)
        :param fasta: Path to a file to write alleles in fasta format (None)
        """
        pass

    @abstractmethod
    def create(self):
        """
        Initialiazes the database for a specific typing method:
            FimH
            Spa
            Clermont
        
        Uses a scheme created temporarily that is specific to the typing, and a multifasta file
        containing target alleles for the method.
        """
        pass


    @abstractmethod
    def multi_search(self, genomes, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        """
        Performed batch search analysis of list of genomes
        
        :param genomes: List of path to the fasta genomes
        :param identity: Minimum identity treshold (0.9)
        :param coverage: Minimum coverage threshold (0.9)
        :param fasta: Path to a file to write alleles in fasta format (None)
        :param output: Write result on this output (stdout)
        """
        pass

    
    def close(self):
        self._database.close()

    def check_input(self, identity, coverage):
        """
        Verify input identity and coverage to be between 0 to 1
        """
        if identity < 0 or identity > 1:
            raise exceptions.BadIdentityRange('Identity must be in range [0-1]')
        logging.info("Use %s minimum identity for the research", str(identity))
        if coverage <0 or coverage > 1:
            raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')
        logging.info("Use %s minimum coverage for the research", str(coverage))

    def write_fasta_allele(self, genome, fasta, psl):
        """
        Export allele in fasta output for a list of psl results
        """
        logging.info("Export allele result to fasta output %s", os.path.basename(fasta.name))
        seqs = read_genome(genome)
        genome_name = file_name(genome)
        for gene,aligns in psl.items():
            for al in aligns:
                seq = seqs.get(al.chro, None)
                if seq is None:
                    raise exceptions.ChromosomeNotFound("Chromosome ID not found " + al.chro)
                al_seq = al.get_sequence(seq)
                fasta.write(">" + genome_name + "|" + gene + "\n" )
                fasta.write(str(al_seq) + "\n")
    

class FimH(PyTyper):
    
    def __init__(self, fi):
        PyTyper.__init__(self, fi, FIM)
        self.typing = FIM

    def search_genome(self, genome):
        print("INSIDE SEARCH GENOME FIM")

    def multi_search(self, genomes):
        print("INSIDE multi_search FIM")
        
    def create(self):
        logging.info("Initialise FIM database")
        with tempfile.TemporaryDirectory() as tmp_dir:
            web.clone_repo(FIM_URL, tmp_dir)
            fas = os.path.join(tmp_dir, "fimH.fsa")
            if os.path.exists(fas) is False:
                raise exceptions.PyMLSTWebError("Problems to retreive fimH fasta file from git repo\n %s", FIM_URL)
            count = 0
            for seq in SeqIO.parse(fas, "fasta"):
                ok = self._database.add_sequence(str(seq.seq), FIM, seq.id)
                if ok:
                    self._database.add_st(seq.id, FIM, seq.id)
                    count += 1
            logging.info("Add %s fimH sequences in the database", str(count))
                

class Spa(PyTyper):
    
    def __init__(self, fi):
        PyTyper.__init__(self, fi, SPA)
        self.typing = SPA

    def search_genome(self, genome, identity, coverage, fasta):
        genome_name = file_name(genome)
        logging.info("Search %s typing for %s genome", self.typing, genome_name)
        result = TypingResult(genome_name, self.typing)
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=True) as tmpout:
            with open(config.get_data('pytyper/spa.fna'),'r') as spafna:
                try:
                    res_all = blat.run_blat(genome, spafna, tmpout, identity, \
                                        0.3, maxintron=1000)
                except exceptions.CoreGenomePathNotFound:
                    logging.warning("No SPA found")
                    result.set_notes("No SPA found")
                    return(result)
        ##check coverage
        res = []
        ids = set()
        for al in res_all.get('spa'):
            if al.coverage > coverage:
                res.append(al)
            else:
                ids.add(al.chro)
        if len(res) == 0:
            logging.info("No complete SPA found")
            result.set_notes("No complete SPA found. Check: " + ";".join(ids))
            return(result)
        else:
            logging.info("Found %s SPA genes", str(len(res)))
            
        ## search individual rep for each potential spa region
        seqs = read_genome(genome)
        notes = []
        count = 0
        for al in res:
            count += 1
            spans = []
            reps = {}
            if al.chro not in seqs:
                raise exceptions.ChromosomeNotFound("Chromosome ID not found " + al.chro)
            sequence = str(al.get_sequence(seqs.get(al.chro)))
            for rep,seq in self._database.get_sequences(self.typing):
                for f in re.finditer(seq, sequence, re.IGNORECASE):
                    spans.append(f.span())
                    reps[f.span()] = rep
            spans.sort()
            ##check found repetition
            if len(spans) == 0:
                continue
            ##check repetition are in a row
            check = self.check_repetition(spans)
            if check is not None:
                notes.append(check)
            ##combinaisons of reps => allele
            allele = "-".join([reps.get(span) for span in spans])
            st = self._database.get_st(self.typing, allele)
            if count==1:
                result.set_st(st)
                result.set_allele(allele)
            else:
                notes.append("Found other Spa : "+st+"|"+allele)
        if len(notes)>0:
            result.set_notes("; ".join(notes))
        if fasta:
            logging.debug("Write alleles in fasta output")
            self.write_fasta_allele(genome, fasta, res)
        return(result)

    def multi_search(self, genomes, identity=0.90, coverage=0.90, fasta=None, output=sys.stdout):
        self.check_input(identity, coverage)
        if coverage<0.8 :
            logging.warning("Use identity less than 0.8 may lead to erroneous results")
        header = True
        for genome in genomes:
            result = self.search_genome(genome, identity, coverage, fasta)
            result.write(output, header)
            if header == True:
                header = False

    def create(self):
        logging.info("Initialise SPA database")
        a = web.request(SPA_URL_TYPE)
        count = 0 
        for line in io.StringIO(a.text):
            h = line.rstrip("\n").split(",")
            if len(h) != 2 :
                logging.warning("Bad line when readind SPA type\n %s", line)
            self._database.add_st(h[0], SPA, h[1].strip())
            count += 1
        logging.info("Add %s SPA type in the database", str(count))
        b = web.request(SPA_URL_SEQ)
        count = 0
        for seq in SeqIO.parse(io.StringIO(b.text), "fasta"):
            ok = self._database.add_sequence(str(seq.seq), SPA, \
                                             seq.id.lstrip("r"))
            if ok:
                count += 1
        logging.info("Add %s SPA sequences in the database", str(count))        
        
    def check_repetition(self, spans):
        """Check that repetition are successive"""
        check = True
        value = None
        notes = []
        for span in spans:
            if value is None:
                value = span[1]
                continue
            if value != span[0]:
                logging.debug("Found discontinous spa repeat at %s", str(value))
                notes.append("Discontinous at "+value)
                check=False
            value = span[1]
        if check:
            return(None)
        else:
            return(", ".append(notes))

class Clmt(PyTyper):
    
    def __init__(self, fi):
        PyTyper.__init__(self, fi, CLMT)
        self.typing = CLMT

    def search_genome(self, genome, identity, coverage, fasta):
        genome_name = file_name(genome)
        logging.info("Search %s typing for %s genome", self.typing, genome_name)
        result = TypingResult(genome_name, self.typing)
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=True) as tmpout:
            with open(config.get_data('pytyper/clmt.fna'),'r') as clmtfna:
                try:
                    res = blat.run_blat(genome, clmtfna, tmpout, identity, coverage)
                    logging.info(", ".join(res.keys()))
                except exceptions.CoreGenomePathNotFound:
                    logging.error("Not Escherichia genome")
                    result.set_notes("Not Escherichia genome")
                    return(result)
        ##check trpA/TrpBA and other escherichia
        if 'trpA' not in res and 'trpBA' not in res:
            result.set_st("Unknown")
            result.set_notes("trpA/trpBA not found")
        elif 'chuAalbertii' in res:
            result.set_allele("chuAalbertii")
            result.set_notes("Escherichia albertii")
        elif 'citPfergus' in res:
            result.set_allele("citPfergus")
            result.set_notes("Escherichia fergusonii")
        else:
            ##check combinaison 4 genes
            als = []
            for gene in ['arpA', 'chuA', 'yjaA', 'TspE4.C2']:
                if gene in res:
                    als.append(gene + "|+")
                else:
                    als.append(gene + "|-")
            result.set_allele(", ".join(als))
            st = self._database.get_st(self.typing, ",".join(als))
            if st is None:
                raise exceptions.AlleleSequenceNotFound(",".join(als) +  " Not found in the database")
            ##check others gene for distinction
            if st == 'I|A|C':
                if 'aesI' in res:
                    result.set_notes("aesI'|+")
                    result.set_st("Clade I")
                elif 'trpAgpC' in res:
                    result.set_notes("aesI'|-, trpAgpC|+")
                    result.set_st("C")
                else:
                    result.set_notes("aesI'|-, trpAgpC|-")
                    result.set_st("A")
            elif st == 'D|E':
                if 'fdm' in res:
                    result.set_notes("fdm|+")
                    result.set_st("E")
                else:
                    result.set_notes("fdm|-")
                    result.set_st("D")
            elif st == 'E|I':
                if 'aesI' in res:
                    result.set_notes("aesI|+")
                    result.set_st("Clade I")
                elif 'fdm' in res:
                    result.set_notes("aesI|-, fdm|+")
                    result.set_st("E")
                else:
                    result.set_notes("aesI|-, fdm|-")
                    result.set_st("E or Clade I")
            elif st == 'I|II':
                if 'aesI' in res:
                    result.set_notes("aesI|+")
                    result.set_st("Clade I")
                elif 'aesII' in res:
                    result.set_notes("aesI|-, aesII|+")
                    result.set_st("Clade II")
                else:
                    result.set_notes("aesI|-, aesII|-")         
                    result.set_st("Unknown")
            elif st == 'G|B2':
                if 'ybgD' in res:
                    result.set_notes("ybgD|+")
                    result.set_st("G")
                else:
                    result.set_notes("ybgD|-")
                    result.set_st("B2")
            elif st == 'G|F':
                if 'ybgD' in res:
                    result.set_notes("ybgD|+")
                    result.set_st("G")
                else:
                    result.set_notes("ybgD|-")
                    result.set_st("F")
            elif st == 'H|B2':
                if 'fdm' in res:
                    result.set_notes("fdm|+")
                    result.set_st("H")
                else:
                    result.set_notes("fdm|-")
                    result.set_st("B2")                    
            elif st == 'III|IV|V':
                note = []
                match = False
                if 'chuIII' in res:
                    note.append("chuIII|+")
                    match = True
                else:
                    note.append("chuIII|-")
                if 'chuIV' in res:
                    note.append("chuIV|+")
                    match = True
                else:
                    note.append("chuIV|-")
                if 'chuV' in res:
                    note.append("chuV|+")
                    match = True
                else:
                    note.append("chuV|-")
                if match:
                    result.set_st("Clade III/IV/V")
                else:
                    result.set_st("Unknown")
                result.set_notes(", ".join(note))
            else:
                result.set_st(st)
            if fasta:
                logging.debug("Write alleles in fasta output")
                self.write_fasta_allele(genome, fasta, res)
        return(result)
        
        
    def multi_search(self, genomes, identity=0.90, coverage=0.99, fasta=None, output=sys.stdout):
        self.check_input(identity, coverage)
        if coverage < 0.98:
            logging.warning("Use identity less than 0.98 may lead to erroneous results")
        header = True
        for genome in genomes:
            result = self.search_genome(genome, identity, coverage, fasta)
            result.write(output, header)
            if header == True:
                header = False
            

    def create(self):
        logging.info("Initialise CLERMONT database")
        fi = config.get_data('pytyper/clmt.fna')
        count = 0
        for seq in SeqIO.parse(fi, "fasta"):
            ok = self._database.add_sequence(str(seq.seq), CLMT, seq.id)
            if ok:
                count += 1
        logging.info("Add %s CLERMONT sequences in the database", str(count))
        count = 0
        with open(config.get_data('pytyper/clmt.txt'), 'r') as alfi:
            header = alfi.readline().rstrip("\n").split(",")
            for line in alfi:
                li = line.rstrip("\n").split(",")
                if len(li) != len(header):
                    raise exceptions.BadInputForCreate("Error when reading clmt profiles\n %s", line)
                alleles = ",".join([h+"|"+v for h,v in zip(header[1:], li[1:])])
                self._database.add_st(li[0], CLMT, alleles)
                count += 1
            logging.info("Add %s CLERMONT type in the database", str(count))


    # def check_result(self, genome, res):
    #     res_ok = {}
    #     seqs = read_genome(genome)
    #     coregenes = read_genome(config.get_data('pytyper/clmt.fna'))
    #     for gene_name,al_all in res.items():
    #         core_seq = coregenes.get(gene_name)
    #         for al in al_all:
    #             if al.coverage == 1:
    #                 res_ok.setdefault(gene_name, []).append(al)
    #             else:
    #                 logging.debug("Align incomplet gene : %s", gene_name)
    #                 seq = seqs.get(al.chro, None)
    #                 if seq is None:
    #                     raise exceptions.ChromosomeNotFound("Chromosome ID not found " + al.chro)
    #                 sequence = al.get_aligned_sequence(seq, core_seq)
    #                 logging.debug(str(sequence))
    #                 logging.debug(str(core_seq))
    #                 if len(sequence) != len(core_seq):
    #                     logging.debug("Partial gene found : %s, skip", gene_name)
    #                 else:
    #                     res_ok.setdefault(gene_name, []).append(al)
    #     return(res_ok)
        


# Result formatting for the different typing methods
class TypingResult:
    """ Writes the results of the TYPING research"""
    
    def __init__(self, genome_name, method):
        """
        :param genome_name: Name of the genome retrieved from the path provided by the user
        :param method: Typing method uses for analysis
        """
        self.name = genome_name
        self.typing = method
        self.st = ""
        self.allele = ""
        self.notes = ""
        
    def write(self,output=sys.stdout, header=True):
        """Writes the results in output file"""
        if header:
            output.write("Sample\tMethod\tType\tAllele\tNotes\n")
        towrite = [self.name, self.typing, self.st, self.allele, self.notes]
        output.write("\t".join(towrite) + "\n")

    def set_allele(self, allele):
        """
        :param allele: Allele of the strain recovered in the core genome
        """
        self.allele = allele
        
    def set_st(self, st_val):
        """
        :param st_val: ST values identified for each genome by search
        """
        self.st = st_val
        
    def set_notes(self, notes):
        """
        :param notes: Add notes to search results
        """ 
        self.notes = notes
        
    def __str__(self):
        return(" ".join([self.name, self.typing, self.st, self.allele]))

