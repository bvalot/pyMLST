import logging
import os
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from alembic import command
from alembic.config import Config
from sqlalchemy import create_engine, inspect, select

from pymlst import config
from pymlst.common import flag, exceptions


def records_to_dict(records):
    seq_dict = {}
    for seq in records:
        seq_dict[seq.id] = seq.seq.upper()
    return seq_dict


def read_genome(handle):
    records = SeqIO.parse(handle, 'fasta')
    return records_to_dict(records)


def write_genome(genome_dict, handle):
    for seq_id, seq in genome_dict.items():
        handle.write('> ' + str(seq_id) + '\n'
                     + str(seq) + '\n')


def strip_file(file):
    found = []  
    if file is not None:
        for line in file:
            found.append(line.rstrip('\n'))
    return found


def compar_seqs(seqs):
    count = 0
    for index in range(0, len(seqs[0])):
        seqs_char = {s[index] for s in seqs}
        if '-' in seqs_char:
            seqs_char.remove('-')
        if len(seqs_char) > 1:
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
    log = config.get_logging_level()
    if log == "DEBUG":
        level = logging.DEBUG
    elif log == "INFO":
        level = logging.INFO
    elif log == "WARNING":
        level = logging.WARNING
    else:
        level = logging.ERROR
    logging.basicConfig(
        level=level,
        format='[%(levelname)s: %(asctime)s] %(message)s')


def clean_kwargs(kwargs):
    """Removes kwargs with None values produced by Click.

    Because of the way the Click library binds
    every arguments and options to kwargs entries,
    when a user doesn't specify an option, its name
    is bound to None in the kwargs dictionary.

    By removing the None entries we can pass the kwargs directly
    to the API core functions without overriding the default values.
    """
    for key, value in kwargs.copy().items():
        if value is None:
            kwargs.pop(key)
    return kwargs

def get_output(kwargs):
    """Extract output from kwargs for extractor
    """
    if 'output' in kwargs:
        out_kwargs = {'output': kwargs['output']}
        kwargs.pop('output')
    else:
        out_kwargs = {}
    return kwargs,out_kwargs


def check_type(conn, mlst_type):
    inspector = inspect(conn)
    tables = inspector.get_table_names()
    if 'mlst_type' not in tables:
        set_type(conn, mlst_type)
        return
    m_t = conn.execute(
        select(flag.mlst_type.c.name)
    ).fetchone()
    if m_t.name != mlst_type:
        raise exceptions.WrongBaseType(
            'The base you are attempting to perform '
            'on belongs to the wrong MLST type')


def set_type(conn, mlst_type):
    flag.metadata.create_all(conn)
    conn.execute(
        flag.mlst_type.insert(),
        name=mlst_type)


def get_updated_engine(path, module):
    env_path = config.get_data(os.path.join('alembic', module))
    alembic_cfg = Config()
    alembic_cfg.set_main_option('script_location', env_path)
    logging.getLogger('alembic').setLevel(logging.CRITICAL)

    if path is None:
        engine = create_engine('sqlite://')  # creates a :memory: database
    else:
        engine = create_engine('sqlite:///' + os.path.abspath(path))

    with engine.begin() as conn:
        check_type(conn, module)
        alembic_cfg.attributes['connection'] = conn
        command.upgrade(alembic_cfg, 'head')

    return engine

def clean_geneid(geneid):
    """Remove '_' on geneid to be compatible with kma search"""
    return(geneid.replace("_", "-"))
