import logging
import os
import gzip
from pathlib import Path
from io import TextIOWrapper

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from alembic import command
from alembic.config import Config
from sqlalchemy import create_engine, inspect, select

from pymlst import config
from pymlst.common import flag, exceptions



def records_to_dict(records):
    """Convert SeqIO records to dictionary and upper seq."""
    seq_dict = {}
    for seq in records:
        seq_dict[seq.id] = seq.seq.upper()
    return seq_dict


def is_fasta_gzip(filepath):
    """Check if Path is a gzip file"""
    if filepath.suffix == '.gz':
        return True
    elif filepath.suffix == '.fasta' or filepath.suffix == '.fna':
        return False
    else:
        with filepath.open(mode='rb') as handle:
            if handle.read(2) == b'\x1f\x8b':
                return True
            else:
                return False
    

def read_genome(filepath):
    """
    Read sequence from an Path object
    
    :param filepath: An file Path object (fasta.(.gz))
    :return: Dictionary with sequence IDs as keys and sequences as values
    """
    if is_fasta_gzip(filepath):
        with gzip.open(filepath, 'rt') as handle:
            records = SeqIO.parse(handle, 'fasta')
            return records_to_dict(records)
    else:
        with filepath.open('rt') as handle:
            records = SeqIO.parse(handle, 'fasta')
            return records_to_dict(records)
    
def handle_fasta_file(filepath):
    """
    Open handle from an Path object
    
    :param filepath: An file Path object (fasta.(.gz))
    """
    if is_fasta_gzip(filepath):
        return gzip.open(filepath, 'rt')
    else:
        return filepath.open('rt')


def write_genome(genome_dict, handle):
    """Write a genome dictionary to a file handle in FASTA format."""
    for seq_id, seq in genome_dict.items():
        handle.write('>' + str(seq_id) + '\n' + str(seq) + '\n')


def strip_file(handle):
    """Read lines from a file handle and strip newline characters."""
    found = []  
    if handle is not None:
        for line in handle:
            found.append(line.rstrip('\n'))
    return found


def compar_seqs(seqs):
    """Compare sequences and count polymorphic sites."""
    count = 0
    for index in range(0, len(seqs[0])):
        seqs_char = {s[index] for s in seqs}
        if '-' in seqs_char:
            seqs_char.remove('-')
        if len(seqs_char) > 1:
            count += 1
    return count


def write_count(count, texte):
    """Write text to a count file if it exists."""
    if count:
        count.write(texte)


def validate_sequence(sequence):
    """Validate that a sequence is a valid coding sequence and doesnt contains ambigus bases"""
    if contains_ambigus(sequence):
        return False
    try:
        sequence.translate(cds=True, table=11)
    except TranslationError:
        return False
    else:
        return True


def contains_ambigus(sequence):
    """Validate that sequence doesnt contains ambigues bases"""
    seqtmp = sequence.upper()
    a = seqtmp.count('A') + seqtmp.count('T') + seqtmp.count('G') + seqtmp.count('C')
    if a != len(sequence):
        return True
    else:
        return False


def create_logger():
    """Create and configure logger based on configuration."""
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
    """Removes kwargs with None values produced by Click."""
    for key, value in kwargs.copy().items():
        if value is None:
            kwargs.pop(key)
    return kwargs


def get_output(kwargs):
    """Extract output from kwargs for extractor."""
    if 'output' in kwargs:
        out_kwargs = {'output': kwargs['output']}
        kwargs.pop('output')
    else:
        out_kwargs = {}
    return kwargs, out_kwargs


def check_type(conn, mlst_type):
    """Check if database matches expected MLST type."""
    inspector = inspect(conn)
    tables = inspector.get_table_names()
    if len(tables) == 0:
        return
    elif 'mlst_type' not in tables:
        logging.warning('The base missing mlst_type metadata, continue with %s', mlst_type)
        return
    m_t = conn.execute(
        select(flag.mlst_type.c.name)
    ).fetchone()
    if m_t.name != mlst_type:
        raise exceptions.WrongBaseType(
            'The base you are attempting to perform '
            'on belongs to the wrong MLST type')


def get_updated_engine(path, module):
    """Get a SQLAlchemy engine with database updated to latest schema."""
    env_path = config.get_data(os.path.join('alembic', module))
    alembic_cfg = Config()
    alembic_cfg.set_main_option('script_location', env_path)
    logging.getLogger('alembic').setLevel(logging.CRITICAL)

    if path is None:
        engine = create_engine('sqlite://')
    else:
        engine = create_engine('sqlite:///' + os.path.abspath(path))

    with engine.begin() as conn:
        check_type(conn, module)
        alembic_cfg.attributes['connection'] = conn
        command.upgrade(alembic_cfg, 'head')

    return engine


def clean_geneid(geneid):
    """Remove '_' on geneid to be compatible with kma search."""
    return geneid.replace("_", "-")
