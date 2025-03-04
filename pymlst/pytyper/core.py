""" Core classes and functions to work in alternative typing methods. """

import logging
import os
import sys
import re
import io
import tempfile

# For type declarations in code
from typing import Union, TypeAlias
File: TypeAlias = Union[str, os.PathLike]


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
from pymlst.common.utils import create_logger

# Creates generator function for pytyper objects (and add context manager decorator) 
@contextmanager
def open_typer(method: str):
    """
    :param method: Defines typing method to apply. Possible values :
        1- fim
        2- spa
        3- clmt

    :yields: A :class: '~pymlst.pytyper.core.pyTyper' object.
    """
    utils.create_logger()
    
    # absolute path to automatically initialized pyTyper database
    fi: File = config.get_data('pytyper.db')

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

    def __init__(self, fi: File) -> None:
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
    def connection(self) -> None:
        return self.__connection

    def check_db(self, method):
        """Checks if database contains specified typing method."""
        count = self.connection.execute(
                select([model.typerSeq.c.id])
                .where(model.typerSeq.c.typing == method)
                ).fetchall()

        return(len(count) > 1)


    def add_sequence(self, sequence, method, allele) -> bool:
        try:
            self.connection.execute(
                model.typerSeq.insert(),
                sequence=sequence, typing=method, allele=allele
                    )
        except IntegrityError:
            logging.warning(f"Duplicate sequence for allele {allele}.")
            return(False)
        return(True)

        
    def add_st(self, st, method, allele) -> None:
        self.connection.execute(
            model.typerSt.insert(),
            st=st, typing=method, allele=allele
        )

    def get_st(self, method, allele) -> str:
        """Gets all the STs of a gene/allele pair."""
        typer_st = self.connection.execute(
            select([model.typerSt.c.st])
            .where(and_(
                model.typerSt.c.typing == method,
                model.typerSt.c.allele == allele))
        ).fetchone()
        return(typer_st[0])
        
    def close(self) -> None:
        self.__connection.close()


# Generic class for pyTyper
class PyTyper(ABC):
    """ Primary class for all pyTyper objects listed on method
    """
    def __init__(self, fi : File, typing : str) -> None:
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

    
    def close(self) -> None:
        self._database.close()

    def check_input(self, identity, coverage):
        if identity < 0 or identity > 1:
            raise exceptions.BadIdentityRange('Identity must be in range [0-1]')
        if coverage <0 or coverage > 1:
            raise exceptions.BadCoverageRange('Coverage must be in range [0-1]')
        

class FimH(PyTyper):
    
    def __init__(self, fi: File) -> None:
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

    def search_genome(self, genome):
        print("INSIDE SEARCH GENOME SPA")

    def multi_search(self, genomes):
        print("INSIDE multi_search SPA")

    def create(self):
        logging.info("Initialise SPA database")
        a = web.request(SPA_URL_TYPE)
        count = 0 
        for line in io.StringIO(a.text):
            h = line.rstrip("\n").split(",")
            if len(h) != 2 :
                logging.warning("Bad line when readind SPA type\n %s", line)
            self._database.add_st(h[0], SPA, h[1])
            count += 1
        logging.info("Add %s SPA type in the database", str(count))
        b = web.request(SPA_URL_SEQ)
        count = 0
        for seq in SeqIO.parse(io.StringIO(b.text), "fasta"):
            ok = self._database.add_sequence(str(seq.seq), SPA, seq.id)
            if ok:
                count += 1
        logging.info("Add %s SPA sequences in the database", str(count))        
        


class Clmt(PyTyper):
    
    def __init__(self, fi):
        PyTyper.__init__(self, fi, CLMT)
        self.typing = CLMT

    def search_genome(self, genome, identity, coverage, fasta):
        genome_name = os.path.basename(genome.name).split('.')[0]
        logging.info("Search %s typing for %s genome", self.typing, genome_name)
        result = TypingResult(genome_name, self.typing)
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=True) as tmpout:
            with open(config.get_data('pytyper/clmt.fna'),'r') as clmtfna:
                try:
                    res = blat.run_blat(genome, clmtfna, tmpout, identity, coverage)
                    logging.debug(", ".join(res.keys()))
                except exceptions.CoreGenomePathNotFound:
                    logging.error("Not Escherichia genome")
                    result.set_notes("Not Escherichia genome")
                    return(result)
        ##check trpA/TrpBA and other escherichia
        if 'trpA' not in res and 'trpBA' not in res:
            result.set_notes("Not Escherichia genome")
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
                if 'Cl.I' in res:
                    result.set_notes("Cl.I'|+")
                    result.set_st("Clade I")
                elif 'Gp.C' in res:
                    result.set_notes("Cl.I'|-, Gp.C|+")
                    result.set_st("C")
                else:
                    result.set_notes("Cl.I'|-, Gp.C|-")
                    result.set_st("A")
            elif st == 'D|E':
                if 'Gp.E' in res:
                    result.set_notes("Gp.E|+")
                    result.set_st("E")
                else:
                    result.set_notes("Gp.E|-")
                    result.set_st("D")
            elif st == 'E|I':
                if 'Cl.I' in res:
                    result.set_notes("Cl.I|+")
                    result.set_st("Clade I")
                elif 'Gp.E' in res:
                    result.set_notes("Cl.I|-, Gp.E|+")
                    result.set_st("E")
                else:
                    result.set_notes("Cl.I|-, Gp.E|-")
                    result.set_st("Unknown")
            elif st == 'I|II':
                if 'Cl.I' in res:
                    result.set_notes("Cl.I|+")
                    result.set_st("Clade I")
                elif 'Cl.II' in res:
                    result.set_notes("Cl.I|-, Cl.II|+")
                    result.set_st("Clade II")
                else:
                    result.set_notes("Cl.I|-, Cl.II|-")         
                    result.set_st("Unknown")
            elif st == 'III|IV|V':
                if 'Cl.III' in res:
                    result.set_notes("Cl.III|+")
                    result.set_st("Clade III")
                elif 'Cl.IV' in res:
                    result.set_notes("Cl.III|-, Cl.IV|+")
                    result.set_st("Clade IV")
                elif 'Cl.V' in res:
                    result.set_notes("Cl.III|-, Cl.IV|-, cl.V|+")
                    result.set_st("Clade V")
                else:
                    result.set_notes("Cl.III|-, Cl.IV|-, cl.V|-")                    
                    result.set_st("Unknown")
            else:
                result.set_st(st)
            if fasta:
                logging.debug("Write alleles in fasta output")
        return(result)
        
        
    def multi_search(self, genomes, identity=0.9, coverage=0.9, fasta=None, output=sys.stdout):
        self.check_input(identity, coverage)
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

