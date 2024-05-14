""" Core classes and functions to work in alternative typing methods. """

import logging
import os
import sys
import re

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
from pymlst.common import blat, kma, utils, exceptions
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
    file: File = os.path.abspath('/home/aimy/Documents/Stage/pyMLST/pymlst/data/alembic/pytyper/idle.db')  # absolute path to automatically initialized pyTyper database

    if method == 'fim':
        typer = FimH(file)
    elif method == 'spa':
        typer = Spa(file)
    elif method == 'clmt':
        typer = Clmt(file)
    else:
        raise exceptions.WrongBaseType(f'Method {method} not implemented.') 
    
    try:
        yield typer
    finally:
        typer.close()






class DatabaseTyper:


    def __init__(self, file: File) -> None:
        """
        :param file: The path to the database file to work with.
        """
        engine = utils.get_updated_engine(file, 'pytyper')
        self.__connection = engine.connect()

        

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
                select([model.typer.c.id])
                .where(model.typer.c.typing == method)
                ).fetchall()

        return(len(count) > 1)


    def add_sequence(self, sequence, method, allele) -> None:
        try:
            self.connection.execute(
                model.typer.insert(),
                sequence=sequence, typing=method, allele=allele
                    )
        except IntegrityError:
            logging.warning(f"Duplicate sequence for allele {allele}.")

    def close(self) -> None:
        self.__connection.close()


# Generic class for pyTyper
class PyTyper(ABC):
    """ Primary class for all pyTyper objects:
        
        - fimH
        - spa
        - clmt

    """
    def __init__(self, file: File) -> None:
        """ 
        :param file: 
        :param ref:
        """
        utils.create_logger()
        self.__database = DatabaseTyper(file)
        self.__file = file

# All abstract methods are declared with empty body, Python interpreter automatically links these methods inside subclasses

    @abstractmethod
    def search_genome(self, genome):
        """
        Abstract method for searching alleles against a genome.
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
    def multi_search(self, genomes):
        pass

    def close(self) -> None:
        self.__database.close()


class FimH(PyTyper):
    
    def __init__(self, file: File) -> None:
        super().__init__(file)
        self.typing = 'fim'
        self.__database = super(PyTyper).__database

    def search_genome(self, genome):
        print("INSIDE SEARCH GENOME FIM")

    def multi_search(self, genomes):
        print(f"{dir(super())}")
        if not self.__database.check_db(self.typing):
            print("EMPTY DATABASE")
        


    def create(self):
         
        self.__tmp = tmp_scheme
        self.__alleles = alleles
        pass
        



class Spa(PyTyper):
    
    def __init__(self, file):
        super().__init__(file)
        self.typing = 'spa'

    def search_genome(self, genome):
        print("INSIDE SEARCH GENOME SPA")


class Clmt(PyTyper):
    
    def __init__(self, file):
        super().__init__(file)
        self.typing = 'clmt'

    def search_genome(self, genome):
        print("INSIDE SEARCH GENOME CLERMONT")


# Result formatting for the different typing methods
class Result:
    pass





